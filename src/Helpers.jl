
# For reference, see https://github.com/frame-consulting/DiffFusion.jl/tree/main/src/utils

"""
    _d_1(strike, forward, nu)

Calculate the Black formula d_1 term.

We implement this once in order to fine-tune corner cases later.
"""
_d_1(strike, forward, nu) = log.(forward ./ strike) ./ nu .+ nu / 2.0


"""
    black_price(strike, forward, nu, call_put)

Calculate Vanilla option price V in Black model.

Argument `strike` (``K``) represents the option strike, `forward` (``F``) is the
forward price of the underlying ``S``, i.e. ``F = E[S(T)]`` in the respective
pricing measure.

Argument nu represents the standard deviation of the forward price. For
annualised lognormal volatility σ, we have ``ν = σ √T``. Finally, `call_put`
is the call (+1) or put (-1) option flag.

Note that we allow negative standard deviation ν. This allows modelling negative
skew for shifted-lognormal models.

Note that we do not test for vanishing volatility. Zero volatility may cause
floating point errors.

We allow broadcasting for arguments.
"""
function black_price(strike, forward, nu, call_put)
    d1 = _d_1(strike, forward, nu)
    d2 = d1 .- nu
    option_value = call_put * (
        forward .* cdf.(Normal(), call_put .* d1) .-
        strike  .* cdf.(Normal(), call_put .* d2)
        )
    intrinsic_value = max.(call_put * (forward .- strike), 0.0)
    return max.(option_value, intrinsic_value)
end

"""
    black_price(strike, forward, σ, T, call_put)

Calculate Vanilla option price ``V`` in Black model with volatility parameter.
"""
black_price(strike, forward, σ, T, call_put) = black_price(strike, forward, σ * sqrt(T), call_put)

"""
    bachelier_price(strike, forward, nu, call_put)

Calculate Vanilla option price ``V`` in Bachelier model.

Argument `strike` (``K``) represents the option strike, `forward` (``F``) is the
forward price of the underlying ``S``, i.e. ``F = E[S(T)]`` in the respective
pricing measure.

Argument `nu` represents the standard deviation of the forward price. For
annualised normal volatility σ, we have `ν = σ √T`. Finally, `call_put`
is the call (+1) or put (-1) option flag.

Note that we do not test for vanishing volatility. Zero volatility may cause
floating point errors.
"""
function bachelier_price(strike, forward, nu, call_put)
    u = call_put .* (forward .- strike)
    h = u ./ nu
    option_value = u .* cdf.(Normal(), h) .+ nu .* pdf.(Normal(), h)
    intrinsic_value = max.(u, 0.0)
    return max.(option_value, intrinsic_value)
end

"""
    bachelier_price(strike, forward, σ, T, call_put)

Calculate Vanilla option price ``V`` in Bachelier model with volatility parameter.
"""
bachelier_price(strike, forward, σ, T, call_put) = bachelier_price(strike, forward, σ * sqrt(T), call_put)


"""
    bachelier_vega(strike, forward, nu)

Calculate Vanilla option Vega in Bachelier model.

Here, Vega is calculated as ``dV / d ν``.
"""
function bachelier_vega(strike, forward, nu)
    h = (forward .- strike) ./ nu
    return (nu .> 0.0) .* pdf.(Normal(), h)
end

"""
    bachelier_vega(strike, forward, σ, T)

Calculate Vanilla option Vega in Bachelier model with volatility parameter.

Here, Vega is calculated as ``dV / dσ``.
"""
bachelier_vega(strike, forward, σ, T) = bachelier_vega(strike, forward, σ * sqrt(T)) * sqrt(T)


"""
    bachelier_implied_stdev(price, strike, forward, call_put, min_max = (1.0e-4, 6.0e-2))

Calculate the implied normal standard deviation ν from a Bachelier model price.
"""
function bachelier_implied_stdev_numeric(price, strike, forward, call_put, min_max = (1.0e-4, 6.0e-2))
    f(nu) = bachelier_price(strike, forward, nu, call_put) - price
    nu = find_zero(f, min_max, Roots.Brent(), xatol=1.0e-8)
    return nu
end

"""
    bachelier_implied_volatility(price, strike, forward, T, call_put, min_max = (0.0001, 0.02))

Calculate the implied normal volatility σ from a Bachelier model price.
"""
function bachelier_implied_volatility_numeric(price, strike, forward, T, call_put, min_max = (0.0001, 0.02))
    return bachelier_implied_stdev_numeric(price, strike, forward, call_put, min_max .* sqrt(T)) / sqrt(T)
end

