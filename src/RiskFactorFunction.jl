
"""
    initial_values(idx, v0, ds, dv, dw, r_extrapolation)

Derive the model parameters for a segment specified by the index `idx`. 
"""
function initial_values(idx, v0, ds, dv, dw, r_extrapolation)
    if idx == 0
        dw0 = 0.0
        ds0 = 0.0
        dv_ = 0.0
        r_ = dv[begin] / ds[begin]
    else
        dw0 = dw[idx]
        ds0 = ds[idx]
        dv_ = dv[idx]
        if idx < length(dv)
            r_ = (dv[idx+1] - dv[idx]) / (ds[idx+1] - ds[idx])
        else
            r_ = r_extrapolation
        end
    end
    return (dw0, ds0, v0 + dv_, r_)
end


"""
    risk_factor(dw_, v0, ds, dv, dw, r_extrapolation)

Calculate the risk factor S from a given auxiliary parameter w.

This is an internal method. Auxiliary parameter and resulting risk factor are
represented relative to ATM.
"""
function risk_factor(dw_, v0, ds, dv, dw, r_extrapolation)
    @assert dw_ ≥ 0.0
    @assert length(ds) == length(dv)
    @assert length(ds) == length(dw)
    idx = searchsortedlast(dw, dw_)
    (dw0, ds0, v_, r_) = initial_values(idx, v0, ds, dv, dw, r_extrapolation)
    if abs(r_) < r_ε
        step = dw_ - dw0
    else
        step = (exp(r_*(dw_-dw0)) - 1.0) / r_
    end
    ds_ = ds0 + v_ * step
    return ds_    
end


"""
    brownian_factor(ds_, v0, ds, dv, dw, r_extrapolation)

Calculate the auxiliary parameter w from a given risk factor S.

This is an internal method. Risk factor and resulting auxiliary parameter are
represented relative to ATM.
"""
function brownian_factor(ds_, v0, ds, dv, dw, r_extrapolation)
    @assert ds_ ≥ 0.0
    @assert length(ds) == length(dv)
    @assert length(ds) == length(dw)
    idx = searchsortedlast(ds, ds_)
    (dw0, ds0, v_, r_) = initial_values(idx, v0, ds, dv, dw, r_extrapolation)
    if abs(r_) < r_ε
        step = 1.0/v_ * (ds_ - ds0)
    else
        tmp = 1.0 + r_/v_*(ds_ - ds0)
        if tmp ≤ 0.0
            step = Inf
        else
            step = log(tmp) / r_
        end
    end
    dw_ = dw0 + step
    return dw_
end


"""
    risk_factor(m::Model, w)

Calculate the risk factor S from a given auxiliary parameter w.
"""
function risk_factor(m::Model, w)
    if w ≥ m.w0
        dw = w - m.w0
        ds = risk_factor(dw, m.v0, m.dsu, m.dvu, m.dwu, m.rexu)
        return m.s0 + ds
    else
        dw = m.w0 - w
        ds = risk_factor(dw, m.v0, m.dsl, m.dvl, m.dwl, m.rexl)
        return m.s0 - ds
    end
end


"""
    brownian_factor(m::Model, s)

Calculate the auxiliary parameter w from a given risk factor S.
"""
function brownian_factor(m::Model, s)
    if s ≥ m.s0
        ds = s - m.s0
        dw = brownian_factor(ds, m.v0, m.dsu, m.dvu, m.dwu, m.rexu)
        return m.w0 + dw
    else
        ds = m.s0 - s
        dw = brownian_factor(ds, m.v0, m.dsl, m.dvl, m.dwl, m.rexl)
        return m.w0 - dw
    end
end


"""
    local_volatility(ds_, v0, ds, dv, dw, r_extrapolation)

Calculate interpolated volatility parameters.
"""
function local_volatility(ds_, v0, ds, dv, dw, r_extrapolation)
    @assert ds_ ≥ 0.0
    @assert length(ds) == length(dv)
    @assert length(ds) == length(dw)
    idx = searchsortedlast(ds, ds_)
    (dw0, ds0, v_, r_) = initial_values(idx, v0, ds, dv, dw, r_extrapolation)
    vol = v_ + r_ * (ds_ - ds0)
    return vol
end


"""
    local_volatility(m::Model, s)

Calculate interpolated volatility parameters.
"""
function local_volatility(m::Model, s)
    if s ≥ m.s0
        ds = s - m.s0
        lv = local_volatility(ds, m.v0, m.dsu, m.dvu, m.dwu, m.rexu)
        return lv
    else
        ds = m.s0 - s
        lv = local_volatility(ds, m.v0, m.dsl, m.dvl, m.dwl, m.rexl)
        return lv
    end
end


"""
    implied_density(m::Model, s)

Calculate model-implied probability density for a given risk factor value.
"""
function implied_density(m::Model, s)
    w = brownian_factor(m, s)
    v = local_volatility(m, s)
    d = pdf.(Normal(), w) / v
    return d
end
