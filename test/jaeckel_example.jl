
using PiecewiseVanillaModel
using Printf
using Test
pvm = PiecewiseVanillaModel

@testset "Jaeckel calibration example" begin

    data = [
        # moneyness        Black vol (Case I)   Black vol (Case II)
        0.035123777453185  0.642412798191439    0.649712512502887
        0.049095433048156  0.621682849924325    0.629372247414191
        0.068624781300891  0.590577891369241    0.598339248024188
        0.095922580089594  0.553137221952525    0.560748840467284
        0.134078990076508  0.511398042127817    0.518685454812697
        0.18741338653678   0.466699250819768    0.473512707134552
        0.261963320525776  0.420225808661573    0.426434688827871
        0.366167980681693  0.373296313420122    0.378806875802102
        0.511823524787378  0.327557513727855    0.332366264644264
        0.715418426368358  0.285106482185545    0.289407658380454
        1                  0.249328882881654    0.253751752243855
        1.39778339939642   0.228967051575314    0.235378088110653
        1.95379843162821   0.220857187809035    0.235343538571543
        2.73098701349666   0.218762825294675    0.260395028879884
        3.81732831143284   0.218742183617652    0.31735041252779
        5.33579814376678   0.218432406892364    0.368205175099723
        7.45829006788743   0.217198426268117    0.417582432865276
        10.4250740447762   0.21573928902421     0.46323707706565
        14.5719954372667   0.214619929462215    0.504386489988866
        20.3684933182917   0.2141074555437      0.539752566560924
        28.4707418310251   0.21457985392644     0.566370957381163
    ]
    T = 5.0722;

    strikes = data[:,1]
    normal_vols_1 = [
        pvm.normal_volatility(black_vol, strike, 1.0, T)
        for (strike, black_vol) in zip(data[:,1], data[:,2])
    ]
    normal_vols_2 = [
        pvm.normal_volatility(black_vol, strike, 1.0, T)
        for (strike, black_vol) in zip(data[:,1], data[:,3])
    ]

    function test_calibration(strikes, normal_vols, T, rexl, rexu, r0l, r0u, α)
        atm_strike = strikes[11]
        atm_vol    = normal_vols[11]
        #
        rel_strikes = strikes .- atm_strike
        rel_strikes = vcat(rel_strikes[1:10], rel_strikes[12:21])
        normal_vols = vcat(normal_vols[1:10], normal_vols[12:21])
        #
        dsl = reverse(-rel_strikes[1:10])
        dsu = rel_strikes[11:20]
        #
        lmfit_kwargs = (
            autodiff = :forwarddiff,
            maxIter  = 100,
            show_trace = false,
        )
        #
        (m, res) = pvm.calibrated_model_from_smile(
            atm_strike, atm_vol, T, dsl, dsu,
            rel_strikes, normal_vols;
            rexl = rexl, rexu = rexu, r0l = r0l, r0u = r0u,
            α=α, lmfit_kwargs = lmfit_kwargs,
        )
        #
        model_implied_b76_vols = [
            pvm.lognormal_volatility(m, strike)
            for strike in m.s0 .+ rel_strikes
        ]
        input_b76_vols = [
            pvm.lognormal_volatility(v, strike, m.s0, m.T)
            for (v, strike) in zip(normal_vols, m.s0 .+ rel_strikes)
        ]
        @test res.converged
        fit = (model_implied_b76_vols .- input_b76_vols) # in decimals
        s = "Converged: " * string(res.converged) * ", "
        s = s * "Fit: " * string(maximum(abs.(fit))) * "\n"
        s = s * pvm.model_string(m)
        # println(s)
        return m, fit
    end


    function test_formulas(m, strikes, black_vols, rel_tol)
        indent = "  "
        s = "Option prices:\n"
        # header
        s = s * indent * "        s "
        s = s * "                o0"
        s = s * "                o1" 
        s = s * "                o2"
        s = s * "                o3\n"
        err_string = ""
        for (strike, σ) in zip(strikes, black_vols)
            cp = (strike ≥ m.s0) ? 1 : -1
            o0 = pvm.black_price(strike, m.s0, σ*sqrt(m.T), cp)
            o1 = pvm.call_put_option_analytic(m, strike, cp)
            o2 = pvm.call_put_option_black(m, strike, cp)
            o3 = pvm.call_put_option_recursive(m, strike, cp)
            #
            s = s * (@sprintf "  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f" strike o0 o1 o2 o3) * "\n"
            #
            err2 = o2/o1 - 1.0
            err3 = o3/o1 - 1.0
            err_string = err_string * (@sprintf "  %16.8f  %16.8f  %16.8f" strike err2 err3) * "\n"
            #
            @test abs(err2) < rel_tol
            @test abs(err3) < rel_tol
        end
        # println(s)
        # println(err_string)
    end


    @testset "Example case I" begin
        rexl, rexu, r0l, r0u, α = nothing, 0.2, 0.0, 0.0, 0.0
        m, fit = test_calibration(strikes, normal_vols_1, T, rexl, rexu, r0l, r0u, α)
        @test maximum(abs.(fit)) < 1.7e-6
        test_formulas(m, strikes, data[:,2], 6.5e-4)  # high tolerance due to extreme forwards
    end

    @testset "Example case II" begin
        rexl, rexu, r0l, r0u, α = nothing, 0.0, 0.0, 0.0, 0.0
        m, fit = test_calibration(strikes, normal_vols_2, T, rexl, rexu, r0l, r0u, α)
        @test maximum(abs.(fit)) < 7.5e-3
        α = 5.0e-2
        m, fit = test_calibration(strikes, normal_vols_2, T, rexl, rexu, r0l, r0u, α)
        test_formulas(m, strikes, data[:,3], 2.0e-11)
    end

end