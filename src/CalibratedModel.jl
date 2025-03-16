
const v_ε = 1.0e-8  # ensure positive volatility parameters

"""
    positive_volatilities(v0, dv)

Check for resulting positive volatilities and floor volatilities if necessary.
"""
function positive_volatilities(v0, dv)
    if (v0 .≥ v_ε) && all(v0 .+ dv .≥ v_ε)
        return (v0, dv)
    end
    v0_ = max(v0, v_ε)
    dv_ = max.(v0_ .+ dv, v_ε) .- v0_
    return (v0_, dv_)
end

function _adjusted_model(x, m0)
    w0 = 0.5 * x[1] * m0.T  # w0 = 0.5*r0*T
    v0 = exp(x[2])  # ensure positivity
    (v0, dvl) = positive_volatilities(v0, m0.dvl)
    (v0, dvu) = positive_volatilities(v0, m0.dvu)
    return model(m0.s0, v0, w0, m0.T, m0.dsl, m0.dsu, dvl, dvu, rexl = m0.rexl, rexu = m0.rexu)
end


function _calibrated_model_F(x, m0, σ_atm, atm_vega)
    m = _adjusted_model(x, m0)
    c = call_option(m, m.s0)
    p = put_option(m, m.s0)
    y1 = (c - p) / atm_vega
    y2 = (c + p) / atm_vega - 2.0 * σ_atm
    return [ y1, y2 ]
end

"""
    calibrated_model(
        s_atm,
        σ_atm,
        T,
        dsl,
        dsu,
        dvl,
        dvu;
        rexl  = 0.0,
        rexu  = 0.0,
        σ_tol = 1.0e-10,
        k_max = 10,
        l_max = 5,
        backtracking_factor = 0.5,
        )

Create a model which is auto-calibrated to at-the-money forward and normal volatility.
"""
function calibrated_model(
    s_atm,
    σ_atm,
    T,
    dsl,
    dsu,
    dvl,
    dvu;
    rexl  = 0.0,
    rexu  = 0.0,
    σ_tol = 1.0e-10,
    k_max = 10,
    l_max = 5,
    backtracking_factor = 0.5,
    )
    #
    @assert σ_atm > 0.0
    @assert T > 0.0
    @assert length(dsl) > 0
    @assert length(dsu) > 0
    @assert length(dvl) > 0
    @assert length(dvu) > 0
    @assert length(dvl) == length(dsl)
    @assert length(dvu) == length(dsu)
    @assert dsl[begin] > 0.0
    @assert dsu[begin] > 0.0
    if length(dsl) > 1
        @assert all((dsl[begin+1:end] .- dsl[begin:end-1]) .> 0.0)
    end
    if length(dsu) > 1
        @assert all((dsu[begin+1:end] .- dsu[begin:end-1]) .> 0.0)
    end
    #
    atm_vega = sqrt(T) / sqrt(2*π)
    r0 = (dvu[begin] - dvl[begin]) / (dsu[begin] + dsl[begin])
    w0 = 0.5*r0*T
    v0 = σ_atm
    (v0, dvl) = positive_volatilities(v0, dvl)
    (v0, dvu) = positive_volatilities(v0, dvu)
    m0 = model(s_atm, v0, w0, T, dsl, dsu, dvl, dvu, rexl = rexl, rexu = rexu)
    F(x) = _calibrated_model_F(x, m0, σ_atm, atm_vega)
    x = [ r0, log(v0) ]
    y = F(x)
    for k = 1:k_max
        x0 = x
        y0 = y
        # println(y0 * 1.0e4)
        if maximum(abs.(y0)) ≤ σ_tol
            break
        end
        J0 = ForwardDiff.jacobian(F, x0)
        dx = -(J0 \ y0)  # full Newton step
        x = x0 + dx
        y = F(x)
        for l = 1:l_max  # line search
            if maximum(abs.(y)) < maximum(abs.(y0))
                break  # accept step
            end
            dx *= backtracking_factor  # simple backtracking
            x = x0 + dx
            y = F(x)
        end
    end
    return _adjusted_model(x, m0)
end

"""
    calibrated_model_from_slopes(
        s_atm,
        σ_atm,
        T,
        std_devs_lo,
        std_devs_up,
        slope_lo,
        slope_up,
        rexl,
        rexu,
        )

Create an ATM-calibrated model from relative/normalised reference strike
and volatility slope increments.

This method is used internally for smile calibration.
"""
function calibrated_model_from_slopes(
    s_atm,
    σ_atm,
    T,
    std_devs_lo,
    std_devs_up,
    slope_lo,
    slope_up,
    rexl,
    rexu,
    )
    #
    @assert length(std_devs_lo) == length(slope_lo)
    @assert length(std_devs_up) == length(slope_up)
    dsl = σ_atm * sqrt(T) * std_devs_lo
    dsu = σ_atm * sqrt(T) * std_devs_up
    #
    slope_lo = cumsum(slope_lo)
    slope_up = cumsum(slope_up)
    step_lo = vcat(dsl[begin], dsl[begin+1:end] .- dsl[begin:end-1])
    step_up = vcat(dsu[begin], dsu[begin+1:end] .- dsu[begin:end-1])
    dvl = cumsum(slope_lo .* step_lo)
    dvu = cumsum(slope_up .* step_up)
    # println(dvl)
    # println(dvu)
    return calibrated_model(s_atm, σ_atm, T, dsl, dsu, dvl, dvu, rexl = rexl, rexu = rexu)
end


function _model_from_x(x, size_lo, s_atm, σ_atm, T, std_devs_lo, std_devs_up, rexl, rexu)
    return calibrated_model_from_slopes(s_atm, σ_atm, T, std_devs_lo, std_devs_up, x[begin:size_lo], x[size_lo+1:end], rexl, rexu)
end

function _calibrated_model_from_smile_F(x, size_lo, s_atm, σ_atm, T, std_devs_lo, std_devs_up, rel_strikes, σ_smiles, rexl, rexu, α)
    m = _model_from_x(x, size_lo, s_atm, σ_atm, T, std_devs_lo, std_devs_up, rexl, rexu)
    σ_model = [
        normal_volatility(m, rel_strike + m.s0) for rel_strike in rel_strikes
    ]
    F_vols = (σ_model .- σ_smiles) ./ σ_atm
    if α > 0.0
        F_x = vcat(x[begin+1:size_lo], x[size_lo+2:end])
        # ensure α acts independent of number of inputs/outputs
        F_vols = (1.0/max(1, length(F_vols))) .* F_vols
        F_x = (1.0/max(1, length(F_x))) .* F_x
        # apply convex combination
        F_vols = vcat((1.0-α) .* F_vols, α .* F_x)
    end
    return F_vols
end

"""
    calibrated_model_from_smile(
        s_atm,
        σ_atm,
        T,
        std_devs_lo,
        std_devs_up,
        rel_strikes,
        σ_smiles;
        rexl  = 0.0,
        rexu  = 0.0,
        α = 0.0,
        σ_min = 1.0e-4,
        lmfit_kwargs = (
            autodiff = :forwarddiff,
            maxIter  = 10
        ),
        )

Create a model that is calibrated to a set of strikes and implied normal volatilities.

Input implied normal volatilities are `σ_smiles`. Corresponding strikes are `rel_strikes`.
Strikes are represented as off-set to at-the-money (ATM) forward `s_atm`.

ATM implied normal volatility is `σ_atm`. Time to option expiry is `T`.

Lower and upper reference strikes are represented by `std_devs_lo` and `std_devs_up`.
The reference strikes are expressed in terms of standard deviations from ATM. Transformation
into actual strikes is implemented via the scaling `σ_atm*sqrt(T)`.

Lower and upper smile extrapolation is controlled via the slope parameters `rexlo` and
`rexu`.

The calibration method is controlled via the named tuple `lmfit_kwargs`. The arguments
are passed on to the `LsqFit.lmfit(...)` method.
"""
function calibrated_model_from_smile(
    s_atm,
    σ_atm,
    T,
    std_devs_lo,
    std_devs_up,
    rel_strikes,
    σ_smiles;
    rexl  = 0.0,
    rexu  = 0.0,
    α = 0.0,
    σ_min = 1.0e-4,
    lmfit_kwargs = (
        autodiff = :forwarddiff,
        maxIter  = 10
    ),
    )
    #
    @assert length(rel_strikes) ≥ 2
    @assert length(rel_strikes) == length(σ_smiles)
    @assert σ_atm ≥ σ_min
    @assert all(σ_smiles .≥ σ_min)
    F(x) = _calibrated_model_from_smile_F(
        x,
        length(std_devs_lo),
        s_atm, σ_atm, T,
        std_devs_lo, std_devs_up,
        rel_strikes, σ_smiles,
        rexl, rexu, α,
    )
    x0 = zeros(size(std_devs_lo).+size(std_devs_up))
    y0 = F(x0)
    res = LsqFit.lmfit(
        F,
        x0,
        eltype(y0)[];
        lmfit_kwargs...
    )
        #
    m = _model_from_x(res.param, length(std_devs_lo), s_atm, σ_atm, T, std_devs_lo, std_devs_up, rexl, rexu)
    return (model = m, result = res)
end
