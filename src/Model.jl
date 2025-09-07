
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
    @assert length(dsl) > 0  # we may relax this constraint
    @assert length(dsu) > 0
    @assert length(dsl) == length(dvl)
    @assert length(dsu) == length(dvu)
    for k in 2:length(dsl)
        @assert dsl[k] > dsl[k-1]
    end
    for k in 2:length(dsu)
        @assert dsu[k] > dsu[k-1]
    end
    dwl = brownian_grid(v0, dsl, dvl)
    dwu = brownian_grid(v0, dsu, dvu)
    #
    rexl = _extrapolation_slope(s0, v0, dsl, dvl, rexl, -1.0)
    rexu = _extrapolation_slope(s0, v0, dsu, dvu, rexu, +1.0)
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


"""
    _extrapolation_slope(s0, v0, ds, dv, r_extrapolation::Number, upp_or_low)

Set slope based on explicit value.
"""
function _extrapolation_slope(s0, v0, ds, dv, r_extrapolation::Number, upp_or_low)
    return r_extrapolation
end

"""
    _extrapolation_slope(s0, v0, ds, dv, r_extrapolation::Nothing, upp_or_low)

Set slope based on linear extrapolation.
"""
function _extrapolation_slope(s0, v0, ds, dv, r_extrapolation::Nothing, upp_or_low)
    return r_extrapolation
end

"""
    _extrapolation_slope(s0, v0, ds, dv, r_extrapolation::String, upp_or_low)

Decide slope based on string value.
"""
function _extrapolation_slope(s0, v0, ds, dv, r_extrapolation::String, upp_or_low)
    if uppercase(r_extrapolation) == "FLAT"
        return _extrapolation_slope(s0, v0, ds, dv, 0.0, upp_or_low)
    end
    if uppercase(r_extrapolation) == "NOTHING"
        return _extrapolation_slope(s0, v0, ds, dv, nothing, upp_or_low)
    end
    if uppercase(r_extrapolation) == "LINEAR"
        if length(ds) > 1
            rex = (dv[end] - dv[end-1]) / (ds[end] - ds[end-1])
        else
            rex = dv[end] / ds[end]
        end
        return rex
    end
    if uppercase(r_extrapolation) == "LOGNORMAL"
        @assert upp_or_low in (-1.0, 1.0)
        s = s0 + upp_or_low * ds[end]
        @assert s > 0.0
        v = v0 + dv[end]
        rex = upp_or_low * v / s
        return rex
    end
    return 0.0 # fall-back
end