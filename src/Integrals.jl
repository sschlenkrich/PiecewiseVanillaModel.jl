
"""
    vanilla_option_integral(m::Model, strike, payoff, call_or_put)

Calculate the option price for a generic `payoff` function using
numerical integration via `quadgk(...)`.
"""
function vanilla_option_integral(m::Model, strike, payoff, call_or_put)
    @assert call_or_put in (-1.0, 1.0)
    w0 = brownian_factor(m, strike)
    f1 = 1.0/sqrt(2*π*m.T)
    p(w) = f1 * exp(-w^2/2/m.T)
    f(w) = payoff(risk_factor(m, w)) * p(w)
    if call_or_put > 0
        w1 = max(w0, Z_∞*sqrt(m.T))
    else
        w1 = min(w0, -Z_∞*sqrt(m.T))
    end
    I = quadgk(f, w0, w1)[1]
    return call_or_put * I
end


"""
    cum_dist(m::Model, s)

Calculate the implied cumulative distribution function for a
risk factor level `s`.
"""
function cum_dist(m::Model, s)
    w = brownian_factor(m, s)
    return cdf(Normal(), w / sqrt(m.T))
end


"""
    integral_one(w, s0, v0, w0, T, r)

Calculate anti-derivative for `S(w)`.

This is an internal method used vor call and put option calculation.
"""
function integral_one(w, s0, v0, w0, T, r)
    z = w / sqrt(T)
    if abs(r) < r_ε
        f = s0 - v0 * w0
        ν = v0 * sqrt(T)
        return f * cdf(Normal(), z) - ν * pdf(Normal(), z)
    else
        h = r * sqrt(T)
        λ = v0/r - s0
        f = (s0 + λ) * exp(r * (r*T/2 - w0))
        return f * cdf(Normal(), z - h) - λ * cdf(Normal(), z)
    end
end


"""
    integral_two(w, s0, v0, w0, T, r)

Calculate anti-derivative for `S(w)^2`.

This is an internal method used vor power option calculation.
"""
function integral_two(w, s0, v0, w0, T, r)
    z = w / sqrt(T)
    if abs(r) < r_ε
        f = s0 - v0 * w0
        ν = v0 * sqrt(T)
        return (f^2 + ν^2) * cdf(Normal(), z) - (v0^2*w + 2*v0*f)*sqrt(T) * pdf(Normal(), z)
    else
        h = r * sqrt(T)
        λ = v0/r - s0
        #
        s1 = λ^2 * cdf(Normal(), z)
        s2 = (s0 + λ) * λ * exp(r * (r*T/2 - w0)) * cdf(Normal(), z - h)
        s3 = (s0 + λ)^2 * exp(2*r*(r*T - w0)) * cdf(Normal(), z - 2*h)
        return s1 - 2*s2 + s3
    end
end


const Z_∞ = 10.0  # Φ(Z_∞) = 1


"""
    call_option_analytic(m::Model, strike)

Calculate a call option for a given (out-of-the-money) strike.
"""
function call_option_analytic(m::Model, strike)
    @assert strike ≥ m.s0
    intrinsic_value = max(m.s0 - strike, 0.0)
    dw = brownian_factor(m, strike) - m.w0
    #
    strike_offset = strike * (1.0 - cum_dist(m, strike))
    # extrapolation segment
    w_start = m.w0 + max(dw, m.dwu[end])
    w_end = max(w_start, Z_∞ * sqrt(m.T))
    idx = length(m.dwu)
    (dw0, ds0, v_, r) = initial_values(idx, m.v0, m.dsu, m.dvu, m.dwu, m.rexu)
    I0 = integral_one(w_start, m.s0 + ds0, v_, m.w0 + dw0, m.T, r)
    I1 = integral_one(w_end,   m.s0 + ds0, v_, m.w0 + dw0, m.T, r)
    I  = I1 - I0
    if dw ≥ m.dwu[end]
        return max(I - strike_offset, intrinsic_value)  # ensure non-negative time value
    end
    # full segments
    idx = searchsortedlast(m.dwu, dw)
    @assert idx < length(m.dwu)  # otherwise we should have returned earlier
    for k = idx+1:length(m.dwu)-1
        (dw0, ds0, v_, r) = initial_values(k, m.v0, m.dsu, m.dvu, m.dwu, m.rexu)
        I0 = integral_one(m.w0 + m.dwu[k],   m.s0 + ds0, v_, m.w0 + dw0, m.T, r)
        I1 = integral_one(m.w0 + m.dwu[k+1], m.s0 + ds0, v_, m.w0 + dw0, m.T, r)
        I += I1 - I0
    end
    # broken segment
    (dw0, ds0, v_, r) = initial_values(idx, m.v0, m.dsu, m.dvu, m.dwu, m.rexu)
    I0 = integral_one(m.w0 + dw,           m.s0 + ds0, v_, m.w0 + dw0, m.T, r)
    I1 = integral_one(m.w0 + m.dwu[idx+1], m.s0 + ds0, v_, m.w0 + dw0, m.T, r)
    I += I1 - I0
    return max(I - strike_offset, intrinsic_value)  # ensure non-negative time value
end


"""
    put_option_analytic(m::Model, strike)

Calculate a put option for a given (out-of-the-money) strike.
"""
function put_option_analytic(m::Model, strike)
    @assert strike ≤ m.s0
    intrinsic_value = max(strike - m.s0, 0.0)
    dw = m.w0 - brownian_factor(m, strike)
    #
    strike_offset = strike * cum_dist(m, strike)
    # extrapolation segment
    w_start = m.w0 - max(dw, m.dwl[end])
    w_end = min(w_start, -Z_∞ * sqrt(m.T))
    idx = length(m.dwl)
    (dw0, ds0, v_, r) = initial_values(idx, m.v0, m.dsl, m.dvl, m.dwl, m.rexl)
    r = -r # beware the sign
    I0 = integral_one(w_start, m.s0 - ds0, v_, m.w0 - dw0, m.T, r)
    I1 = integral_one(w_end,   m.s0 - ds0, v_, m.w0 - dw0, m.T, r)
    I  = I0 - I1
    if dw ≥ m.dwl[end]
        return max(strike_offset - I, intrinsic_value) # ensure non-negative time value
    end
    # full segments
    idx = searchsortedlast(m.dwl, dw)
    @assert idx < length(m.dwl)  # otherwise we should have returned earlier
    for k = idx+1:length(m.dwl)-1
        (dw0, ds0, v_, r) = initial_values(k, m.v0, m.dsl, m.dvl, m.dwl, m.rexl)
        r = -r # beware the sign
        I0 = integral_one(m.w0 - m.dwl[k],   m.s0 - ds0, v_, m.w0 - dw0, m.T, r)
        I1 = integral_one(m.w0 - m.dwl[k+1], m.s0 - ds0, v_, m.w0 - dw0, m.T, r)
        I += I0 - I1
    end
    # broken segment
    (dw0, ds0, v_, r) = initial_values(idx, m.v0, m.dsl, m.dvl, m.dwl, m.rexl)
    r = -r # beware the sign
    I0 = integral_one(m.w0 - dw,           m.s0 - ds0, v_, m.w0 - dw0, m.T, r)
    I1 = integral_one(m.w0 - m.dwl[idx+1], m.s0 - ds0, v_, m.w0 - dw0, m.T, r)
    I += I0 - I1
    return max(strike_offset - I, intrinsic_value) # ensure non-negative time value
end


"""
    power_call_option_analytic(m::Model, strike)

Calculate a power call option for a given (out-of-the-money) strike.
"""
function power_call_option_analytic(m::Model, strike)
    @assert strike ≥ m.s0
    intrinsic_value = 0.0  # out-of-the-money
    dw = brownian_factor(m, strike) - m.w0
    #
    strike_offset_1 = 2 * strike * call_option(m::Model, strike)
    strike_offset_2 = strike^2 * (1.0 - cum_dist(m, strike))
    strike_offset   = strike_offset_1 + strike_offset_2
    # extrapolation segment
    w_start = m.w0 + max(dw, m.dwu[end])
    w_end = max(w_start, Z_∞ * sqrt(m.T))
    idx = length(m.dwu)
    (dw0, ds0, v_, r) = initial_values(idx, m.v0, m.dsu, m.dvu, m.dwu, m.rexu)
    I0 = integral_two(w_start, m.s0 + ds0, v_, m.w0 + dw0, m.T, r)
    I1 = integral_two(w_end,   m.s0 + ds0, v_, m.w0 + dw0, m.T, r)
    I  = I1 - I0
    if dw ≥ m.dwu[end]
        return max(I - strike_offset, intrinsic_value)  # ensure non-negative time value
    end
    # full segments
    idx = searchsortedlast(m.dwu, dw)
    @assert idx < length(m.dwu)  # otherwise we should have returned earlier
    for k = idx+1:length(m.dwu)-1
        (dw0, ds0, v_, r) = initial_values(k, m.v0, m.dsu, m.dvu, m.dwu, m.rexu)
        I0 = integral_two(m.w0 + m.dwu[k],   m.s0 + ds0, v_, m.w0 + dw0, m.T, r)
        I1 = integral_two(m.w0 + m.dwu[k+1], m.s0 + ds0, v_, m.w0 + dw0, m.T, r)
        I += I1 - I0
    end
    # broken segment
    (dw0, ds0, v_, r) = initial_values(idx, m.v0, m.dsu, m.dvu, m.dwu, m.rexu)
    I0 = integral_two(m.w0 + dw,           m.s0 + ds0, v_, m.w0 + dw0, m.T, r)
    I1 = integral_two(m.w0 + m.dwu[idx+1], m.s0 + ds0, v_, m.w0 + dw0, m.T, r)
    I += I1 - I0
    return max(I - strike_offset, intrinsic_value)  # ensure non-negative time value
end


"""
    power_put_option_analytic(m::Model, strike)

Calculate a power put option for a given (out-of-the-money) strike.
"""
function power_put_option_analytic(m::Model, strike)
    @assert strike ≤ m.s0
    intrinsic_value = 0.0  # out-of-the-money
    dw = m.w0 - brownian_factor(m, strike)
    #
    strike_offset_1 = strike^2 * cum_dist(m, strike)
    strike_offset_2 = 2 * strike * put_option(m::Model, strike)
    strike_offset   = strike_offset_1 - strike_offset_2
    # extrapolation segment
    w_start = m.w0 - max(dw, m.dwl[end])
    w_end = min(w_start, -Z_∞ * sqrt(m.T))
    idx = length(m.dwl)
    (dw0, ds0, v_, r) = initial_values(idx, m.v0, m.dsl, m.dvl, m.dwl, m.rexl)
    r = -r # beware the sign
    I0 = integral_two(w_start, m.s0 - ds0, v_, m.w0 - dw0, m.T, r)
    I1 = integral_two(w_end,   m.s0 - ds0, v_, m.w0 - dw0, m.T, r)
    I  = I0 - I1
    if dw ≥ m.dwl[end]
        return max(I - strike_offset, intrinsic_value) # ensure non-negative time value
    end
    # full segments
    idx = searchsortedlast(m.dwl, dw)
    @assert idx < length(m.dwl)  # otherwise we should have returned earlier
    for k = idx+1:length(m.dwl)-1
        (dw0, ds0, v_, r) = initial_values(k, m.v0, m.dsl, m.dvl, m.dwl, m.rexl)
        r = -r # beware the sign
        I0 = integral_two(m.w0 - m.dwl[k],   m.s0 - ds0, v_, m.w0 - dw0, m.T, r)
        I1 = integral_two(m.w0 - m.dwl[k+1], m.s0 - ds0, v_, m.w0 - dw0, m.T, r)
        I += I0 - I1
    end
    # broken segment
    (dw0, ds0, v_, r) = initial_values(idx, m.v0, m.dsl, m.dvl, m.dwl, m.rexl)
    r = -r # beware the sign
    I0 = integral_two(m.w0 - dw,           m.s0 - ds0, v_, m.w0 - dw0, m.T, r)
    I1 = integral_two(m.w0 - m.dwl[idx+1], m.s0 - ds0, v_, m.w0 - dw0, m.T, r)
    I += I0 - I1
    return max(I - strike_offset, intrinsic_value) # ensure non-negative time value
end


"""
    variance_analytic(m::Model)

Calculate the model-implied variance of the risk factor.
"""
function variance_analytic(m::Model)
    return power_put_option_analytic(m, m.s0) + power_call_option_analytic(m, m.s0)
end


"""
    call_put_option_analytic(m::Model, strike, call_or_put)

Calculate a call or put option for a given (out-of-the-money) strike.
"""
function call_put_option_analytic(m::Model, strike, call_or_put)
    @assert call_or_put in (-1, 1)
    if call_or_put == 1
        return call_option_analytic(m, strike)
    else
        return put_option_analytic(m, strike)
    end
end