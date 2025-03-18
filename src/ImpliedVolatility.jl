
"""
    lognormal_volatility(m::Model, strike)

Calculate the Black or log-normal implied volatility for a given model.
"""
function lognormal_volatility(m::Model, strike)
    if strike ≥ m.s0
        cp = 1.0
        o = call_option(m, strike)
    else
        cp = -1.0
        o = put_option(m, strike)
    end
    v = black_implied_volatility(o, strike, m.s0, m.T, cp)
    return v
end


"""
    normal_volatility(m::Model, strike)

Calculate the Bachelier or normal implied volatility for a given model.
"""
function normal_volatility(m::Model, strike)
    if strike ≥ m.s0
        cp = 1.0
        o = call_option(m, strike)
    else
        cp = -1.0
        o = put_option(m, strike)
    end
    v = bachelier_implied_volatility(o, strike, m.s0, m.T, cp)
    return v
end


"""
    normal_volatility(b76_volatility, strike, forward, T)

Calculate Bachelier or normal implied volatility from Black or log-normal volatility
"""
function normal_volatility(b76_volatility, strike, forward, T)
    cp = (strike ≥ forward) ? 1.0 : -1.0
    o = black_price(strike, forward, b76_volatility*sqrt(T), cp)
    v = bachelier_implied_volatility(o, strike, forward, T, cp)
    return v
end
