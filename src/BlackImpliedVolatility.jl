
# This file contains methods to calculate implied log-normal volatilities from
# https://github.com/jherekhealy/AQFED.jl.
#
# We include the code here directly to avoid the dependency on the unregistered
# AQFED.jl package.
#
# For code reference, see
#
# https://github.com/jherekhealy/AQFED.jl/blob/master/src/black/iv_solver_common.jl
# https://github.com/jherekhealy/AQFED.jl/blob/master/src/black/iv_solver_householder.jl
#


function black_implied_volatility(price, strike, forward, T, call_put)
    isCall = (call_put==1.0)
    df = 1.0
    ftolrel = 0.0  # see AQFED.jl
    maxEval = 64   # see AQFED.jl
    return impliedVolatilitySRHalley(isCall, price, forward, strike, T, df, ftolrel, maxEval)
end


function normalizePrice(isCall::Bool, price::TP, f::T, strike::T, df::T) where {T,TP}
    c = price / f / df
    ex = f / strike

    if !isCall
        if ex <= 1
            c = c + 1 - 1 / ex # put call parity
        else
            #else duality + put call parity
            c = ex * c
            ex = 1 / ex
        end
    else
        if ex > 1
            # use c(-x0, v)
            c = (f * (c - 1) + strike) / strike # in out duality, c = ex*c + 1 - ex  //not as good numericall
            ex = 1 / ex
        end
    end
    return c, ex
end

const SqrtEpsilon = sqrt(eps())

# f = 1.0, strike = 1/ex, ex=ey = 1.0/strike=f/strike, y = log(f/strike) = log(ey), alpha = price*ex
#polya is 2 / pi, polya_factor = sqrt(pi/8) is aludaat,
function impliedVolatilitySRGuessUndiscountedCall(price::Real, ey::Real, y::Real, polya_factor::Real = 2/pi)::Real
    alpha = price *ey
    r = 2 * alpha - ey + 1
    em2piy = exp(-polya_factor * y)
    A = (ey * em2piy - 1.0 / (ey * em2piy))^2
    r2 = r * r
    B = 4 * (em2piy + 1.0 / em2piy) - 2 / ey * (ey * em2piy + 1.0 / (ey * em2piy)) * ((ey * ey) + 1 - r2)
    beta = 0.0
    if abs(alpha) < SqrtEpsilon
        C = -16 * (1 - 1.0 / ey) * alpha - 16 * (1 - 3 / ey + 1 / (ey * ey)) * alpha * alpha
        beta = C / B
        #    println(" alpha ",alpha, " C ",C," beta ",beta)
    elseif abs(ey - alpha) < SqrtEpsilon
        C = -16 * (1 + 1.0 / ey) * (alpha - ey) - 16 * (1 + 3 / ey + 1 / (ey * ey)) * (alpha - ey) * (alpha - ey)
        if C == 0
            beta = B / A #this is wrong sign but gives order of magnitude
        else
            beta = C / B
        end
        #    println(" alpha ",alpha, " C ",C," beta ",beta)
    elseif abs(y) < SqrtEpsilon
        a2 = alpha * alpha
        a4 = a2 * a2
        C = 16 * ((a2 - a4) + (2 * a4 + 2 * alpha * a2 - a2 - alpha) * y)
        B = 16 * (a2 - y * (alpha + a2))
        beta = C / B
    else
        eym1 = ey - 1
        eyp1 = ey + 1
        C = 1.0 / (ey * ey) * (r2 - eym1 * eym1) * (eyp1 * eyp1 - r2)
        beta = 2 * C / (B + sqrt(B * B + 4 * A * C))
        # println("beta ",beta, " ",C, " ",B, " ",A)
    end
    if beta < 0
        beta = SqrtEpsilon
    end
    gamma = -log(beta) / polya_factor
    if y >= 0
        Asqrty = 0.5 * (1 + sqrt(1 - em2piy * em2piy))
        c0 = (Asqrty - 0.5/ey)
        gmy = max(gamma - y,0)      #machine epsilon issues
        if price <= c0
            return sqrt(gamma + y) - sqrt(gmy)
        end
        return (sqrt(gamma + y) + sqrt(gmy))
    end
    Asqrty = 0.5 * (1 - sqrt(1 - 1.0 / (em2piy * em2piy)))
    c0 = (0.5 - Asqrty/ey)
    gpy = gamma + y
    if gpy < 0
        #machine epsilon issues
        gpy = 0.0
    end
    if price <= c0
        return (-sqrt(gpy) + sqrt(gamma - y))
    end
    return (sqrt(gpy) + sqrt(gamma - y))
end

function impliedVolatilitySRGuess(isCall::Bool, price::Real, f::Real, strike::Real, tte::Real, df::Real)::Real
    c,ex = normalizePrice(isCall, price, f, strike, df)
    scaledVol = impliedVolatilitySRGuessUndiscountedCall(c, ex, log(ex))
    return scaledVol / sqrt(tte)
end

function objectiveHouseholderLog(x::T, ex::T, v::TP, logc::TP) where {T,TP}
    v = abs(v)
    h = x / v
    t = v / 2
    h2 = h^2
    t2 = t^2
    sqrt2 = sqrt(T(2))
    Np = erfcx(-(h + t) / sqrt2)
    Nm = erfcx(-(h - t) / sqrt2)
    eh2t2 = exp(-(h2 + t2) / 2)
    norm = 1 / (2 * sqrt(ex)) * eh2t2
    cEstimate = norm * (Np - Nm)
    logcEstimate = log(cEstimate)
    twopi = 2 * T(pi)
    logvega = (2 /sqrt(twopi)) / (Np - Nm)  #u'/u
    volgaOverVega = (h + t) * (h - t) / v
    logvolgaOverVega = volgaOverVega - logvega  #u'' u - u'2 /u^2
    h2mt2 = h2 - t2
    c3OverVega = (-3 * h2 - t2 + h2mt2 * h2mt2) / (v^2)   #
    logc3overVega = c3OverVega - 3 * logvega * volgaOverVega + 2 * logvega * logvega
    return logcEstimate - logc, (logcEstimate - logc) / logvega, logvolgaOverVega, logc3overVega
end

function impliedVolatilitySRHalley(
    isCall::Bool,
    price::TP,
    f::T,
    strike::T,
    tte::T,
    df::T,
    ftolrel::T,
    maxEval::Int,
    # solver::Householder,  # no need to specialise here
)::TP where {T,TP}
    c, ex = normalizePrice(isCall, price, f, strike, df)
    if c >= 1 / ex || c <= 0
        throw(DomainError(c, string("Price out of range, must be < ", 1 / ex, " and > 0. denormalized price=",price)))
    end
    x = log(ex)
    guess = TP(impliedVolatilitySRGuessUndiscountedCall(c, ex,x))
    b = guess
    xtolrel = 32 * eps(T)
    xtol = T(0)
    ftolrel = max(ftolrel, eps(T))
    ftol = ftolrel * max(1, c)
    logc = log(c)
    fb, fbOverfpb, fp2bOverfpb, fp3bOverfpb = objectiveHouseholderLog(x, ex, b, logc)
    if abs(fb) < ftol
        return b / sqrt(tte)
    end
    for iteration = 1:maxEval
        x0 = b
        lf = (fbOverfpb) * (fp2bOverfpb)
        hn = -fbOverfpb
        num = (1 + fp2bOverfpb * hn / 2)
        denom = (1 + fp2bOverfpb * hn + fp3bOverfpb / 6 * hn * hn)
        x1 = x0 + hn * num / denom
        a = x0
        fb, fbOverfpb, fp2bOverfpb, fp3bOverfpb = objectiveHouseholderLog(x, ex, x1, logc)
        b = x1
        xtol_ = xtol + max(1, abs(b)) * xtolrel
        if abs(b - a) <= xtol_ || abs(fb) < ftol
            break
        end
    end
    return b / sqrt(tte)
end
