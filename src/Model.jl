
const Model = NamedTuple  # an alias

const r_ε = 1.0e-8  # r -> 0, avoid division by zero


"""
    brownian_grid(v0, ds, dv)

Calculate the grid of auxiliary parameters at model creation.
"""
function brownian_grid(v0, ds, dv)
    @assert length(ds) == length(dv)
    dw = zeros(0)
    for k=1:length(ds)
        ds0_ = (k==1) ? 0.0 : ds[k-1]
        dv0_ = (k==1) ? 0.0 : dv[k-1]
        dw0_ = (k==1) ? 0.0 : dw[k-1]
        #
        ds1_ = ds[k]
        dv1_ = dv[k]
        #
        r1 = (dv1_ - dv0_) / (ds1_ - ds0_)
        if abs(r1) < r_ε
            step = 1.0 / (v0 + dv0_)
        else
            step = (log(v0 + dv1_) - log(v0 + dv0_)) / (dv1_ - dv0_)
        end
        dw_k = dw0_ + (ds1_ - ds0_) * step
        dw = vcat(dw, dw_k)
    end
    return dw
end


"""
    model(s0, v0, w0, T, dsl, dsu, dvl, dvu; rexl = 0.0, rexu = 0.0)

Create a model from direct inputs and without calibration.

Model attributes are as follows:

`s0` is the forward risk factor (at-the-money level).

`v0` is the volatility parameter at `s0`.

`w0` is the auxiliary parameter at `s0`.

`T` is the time to option expiry.

`dsl` is a vector of relative reference strikes for lower smile; Actual strikes are `sl = s0 - dsl`.

`dsu` is a vector of relative reference strikes for upper smile; Actual strikes are `su = s0 + dsu`.

`dvl` is a vector of volatility parameter offsets for lower smile; Actual volatility is `vl = v0 + dvl`.

`dvu` is a vector of volatility parameter offsets for upper smile; Actual volatility is `vu = v0 + dvu`.

`rexl` is the outward slope of the lower extrapolation smile; log-normal model uses `rexl < 0`.
If `rexl == nothing` then slope is calculated via linear extrapolation.

`rexu` is the outward slope of the upper extrapolation smile; log-normal model uses `rexl > 0`.
If `rexu == nothing` then slope is calculated via linear extrapolation.
"""
function model(s0, v0, w0, T, dsl, dsu, dvl, dvu; rexl = 0.0, rexu = 0.0)
    @assert v0 > 0.0
    @assert T > 0.0
    @assert length(dsl) == length(dvl)
    @assert length(dsu) == length(dvu)
    # TODO: test monotonicity of S-grids, positivity of v
    dwl = brownian_grid(v0, dsl, dvl)
    dwu = brownian_grid(v0, dsu, dvu)
    #
    if isnothing(rexl)  # low-strike extrapolation
        if length(dsl) > 1
            rexl = (dvl[end] - dvl[end-1]) / (dsl[end] - dsl[end-1])
        else
            rexl = dvl[end] / dsl[end]
        end
    end
    if isnothing(rexu)  # high-strike extrapolation
        if length(dsu) > 1
            rexu = (dvu[end] - dvu[end-1]) / (dsu[end] - dsu[end-1])
        else
            rexu = dvu[end] / dsu[end]
        end
    end
    #
    return (
        s0  = s0,
        v0  = v0,
        w0  = w0,
        T   = T,
        dsl = dsl,
        dsu = dsu,
        dvl = dvl,
        dvu = dvu,
        dwl = dwl,
        dwu = dwu,
        rexl = rexl,
        rexu = rexu,
    )
end
