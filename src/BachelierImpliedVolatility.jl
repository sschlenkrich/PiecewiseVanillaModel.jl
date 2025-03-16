

_Φ_tilde(x) = cdf(Normal(), x) + pdf(Normal(), x) / x

"""
    bachelier_implied_volatility(price, strike, forward, T, call_put)

Calculate the Bachelier or normal implied volatility.

`price` is the forward or undiscounted option price, `strike` is the option strike,
`forward` is the forward asset price or expectation of factor, `T` is time to option
expiry, and `call_put` encodes call (`+1`) or put (`-1`) options.
"""
function bachelier_implied_volatility(price, strike, forward, T, call_put)
    # P. Jäckel, Implied Normal Volatility, 2017
    if isapprox(strike, forward)  # ATM option
        atm_vega = sqrt(T) / sqrt(2*π)
        return price / atm_vega
    end
    time_value = price - max(call_put*(forward-strike), 0.0)
    if isapprox(time_value, 0.0)  # option = intrinsic value
        return 0.0
    end
    @assert	time_value > 0.0  "time_value: " * string(time_value)
    Φ_tilde_star = -abs(time_value) / abs(forward - strike)
    # solve Φ_tilde(x_star) = Φ_tilde_star
    if Φ_tilde_star < -0.001882039271
        g = 1.0 / (Φ_tilde_star - 0.5)
        ξ_bar = 
            (0.032114372355 -
             g * g *
                 (0.016969777977 - g * g * (2.6207332461E-3 - 9.6066952861E-5 * g * g))) /
            (1.0 -
             g * g * (0.6635646938 - g * g * (0.14528712196 - 0.010472855461 * g * g)))
        x_bar = g * (0.3989422804014326 + ξ_bar * g * g)
    else
        h = sqrt(-log(-Φ_tilde_star))
        x_bar =
            (9.4883409779 - h * (9.6320903635 - h * (0.58556997323 + 2.1464093351 * h))) /
            (1.0 - h * (0.65174820867 + h * (1.5120247828 + 6.6437847132E-5 * h)))
    end
    q = (_Φ_tilde(x_bar) - Φ_tilde_star) / pdf(Normal(), x_bar)
    x_star = x_bar +
        3.0 * q * x_bar * x_bar * (2.0 - q * x_bar * (2.0 + x_bar * x_bar)) /
        (6.0 + q * x_bar * (-12.0 + x_bar * (6.0 * q + x_bar *
        (-6.0 + q * x_bar * (3.0 + x_bar * x_bar)))))
    σ = abs((forward - strike)/x_star) / sqrt(T)
    return σ
end
