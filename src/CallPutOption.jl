
@enum VanillaModelType bachelier shiftedlognormal

function model_parameters(m::Model, strike, call_or_put)
    @assert call_or_put in (-1, 1)
    # We only model ATM and OTM options. And for ATM we need an additional flag
    @assert (call_or_put != 1)  || (strike ≥ m.s0)  # call => K ≥ S0
    @assert (call_or_put != -1) || (strike ≤ m.s0)  # put => K ≤ S0
    ds = call_or_put * (strike - m.s0)
    ds = ds + eps()  # make sure we get the outer segment
    if call_or_put == 1
        idx = searchsortedlast(m.dsu, ds)
        (dw0, ds0, v, r) = initial_values(idx, m.v0, m.dsu, m.dvu, m.dwu, m.rexu)        
    else
        idx = searchsortedlast(m.dsl, ds)
        (dw0, ds0, v, r) = initial_values(idx, m.v0, m.dsl, m.dvl, m.dwl, m.rexl)
    end
    w0 = m.w0 + call_or_put * dw0
    s0 = m.s0 + call_or_put * ds0
    r = call_or_put * r  # beware the sign
    if abs(r) < r_ε
        m_type = bachelier
        σ = v
        λ = 0.0
        f = s0 - v * w0
    else
        m_type = shiftedlognormal
        σ = r
        λ = v / r - s0
        f = (s0 + λ) * exp(r*(r*m.T/2.0 - w0)) - λ
    end
    return (m_type=m_type, σ=σ, λ=λ, f=f, idx=idx)
end

function option_value(m_type, strike, ν, λ, f, call_or_put)
    @assert call_or_put in (-1, 1)
    if m_type == shiftedlognormal
        o = black_price(strike + λ, f + λ, ν, call_or_put)
    else
        o = bachelier_price(strike, f, ν, call_or_put)
    end
    return o
end

function call_put_option_recursive(m::Model, strike, call_or_put)
    @assert call_or_put in (-1, 1)
    # We only model ATM and OTM options. And for ATM we need an additional flag
    @assert (call_or_put != 1)  || (strike ≥ m.s0)  # call => K ≥ S0
    @assert (call_or_put != -1) || (strike ≤ m.s0)  # put => K ≤ S0
    p0 = model_parameters(m, strike, call_or_put)
    V0 = option_value(p0.m_type, strike, p0.σ * sqrt(m.T), p0.λ, p0.f, call_or_put)
    if call_or_put == 1
        ds = m.dsu
    else
        ds = m.dsl
    end
    if p0.idx == length(ds)
        return V0
    end
    @assert p0.idx < length(ds)  # otherwise we returned
    s1 = m.s0 + call_or_put * ds[p0.idx+1]
    V1 = option_value(p0.m_type, s1, p0.σ * sqrt(m.T), p0.λ, p0.f, call_or_put)
    return (V0 - V1) + call_put_option_recursive(m::Model, s1, call_or_put)
end


function call_put_option_black(m::Model, strike, call_or_put)
    @assert call_or_put in (-1, 1)
    # We only model ATM and OTM options. And for ATM we need an additional flag
    @assert (call_or_put != 1)  || (strike ≥ m.s0)  # call => K ≥ S0
    @assert (call_or_put != -1) || (strike ≤ m.s0)  # put => K ≤ S0
    p0 = model_parameters(m, strike, call_or_put)
    V = option_value(p0.m_type, strike, p0.σ * sqrt(m.T), p0.λ, p0.f, call_or_put)
    ds = (call_or_put == 1) ? (m.dsu) : (m.dsl)
    for idx in (p0.idx+1):length(ds)
        s1 = m.s0 + call_or_put * ds[idx]
        V0 = option_value(p0.m_type, s1, p0.σ * sqrt(m.T), p0.λ, p0.f, call_or_put)
        p1 = model_parameters(m, s1, call_or_put)
        V1 = option_value(p1.m_type, s1, p1.σ * sqrt(m.T), p1.λ, p1.f, call_or_put)
        V = V + (V1 - V0)
        p0 = p1
    end
    intrinsic_value = max(call_or_put*(m.s0 - strike), 0.0)
    return max(V, intrinsic_value)
end
