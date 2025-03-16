
using PiecewiseVanillaModel
using Distributions
using Test
pvm = PiecewiseVanillaModel

@testset "Calibrated model" begin

    function test_implied_vols(m, strikes)
        for strike in strikes
            cp = (strike≥m.s0) ? 1.0 : -1.0
            option_price = (strike≥m.s0) ? pvm.call_option : pvm.put_option
            v = pvm.normal_volatility(m, strike)
            o1 = pvm.bachelier_price(strike, m.s0, v, m.T, cp)
            o2 = option_price(m, strike)
            @test v ≥ 0.0
            @test isapprox(o1, o2)
        end
    end


    @testset "Upward skew model" begin
        s_atm = 0.03
        σ_atm = 0.01
        r = 0.2
        dsl = [ 0.01, 0.02, 0.04, 0.05 ]
        dsu = [ 0.01, 0.02, 0.03 ]
        dvl = -r .* dsl
        dvu = [ 0.0020, 0.0020, 0.0030 ]
        for T in [ 1.0/12, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0 ]
            m = pvm.calibrated_model(s_atm, σ_atm, T, dsl, dsu, dvl, dvu, rexl = -r, rexu = 0.0)
            c = pvm.call_option(m, m.s0)
            p = pvm.put_option(m, m.s0)
            #
            atm_vega = sqrt(T) * pdf(Normal(), 0.0)
            @test isapprox(c/atm_vega, p/atm_vega, atol=1.0e-10 )
            @test isapprox((c+p)/atm_vega, 2*σ_atm, atol=1.0e-10 )
            test_implied_vols(m, m.s0 .- m.dsl)
            test_implied_vols(m, m.s0 .+ m.dsu)
            #
            m = pvm.calibrated_model(s_atm, σ_atm, T, dsl, dsu, dvl, dvu, rexl = nothing, rexu = nothing)
            c = pvm.call_option(m, m.s0)
            p = pvm.put_option(m, m.s0)
            #
            atm_vega = sqrt(T) * pdf(Normal(), 0.0)
            @test isapprox(c/atm_vega, p/atm_vega, atol=1.0e-10 )
            @test isapprox((c+p)/atm_vega, 2*σ_atm, atol=1.0e-10 )
            test_implied_vols(m, m.s0 .- m.dsl)
            test_implied_vols(m, m.s0 .+ m.dsu)
        end
    end

    @testset "Symmetric down-ward skew model" begin
        s_atm = 0.03
        σ_atm = 0.01
        r = 0.2
        dsu = [ 0.01, 0.02, 0.04, 0.05 ]
        dsl = [ 0.01, 0.02, 0.03 ]
        dvu = -r .* dsu
        dvl = [ 0.0020, 0.0020, 0.0030 ]
        for T in [ 1.0/12, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0 ]
            m = pvm.calibrated_model(s_atm, σ_atm, T, dsl, dsu, dvl, dvu, rexl = r, rexu = -r)
            c = pvm.call_option(m, m.s0)
            p = pvm.put_option(m, m.s0)
            #
            atm_vega = sqrt(T) * pdf(Normal(), 0.0)
            @test isapprox(c/atm_vega, p/atm_vega, atol=1.0e-10 )
            @test isapprox((c+p)/atm_vega, 2*σ_atm, atol=1.0e-10 )
            test_implied_vols(m, m.s0 .- m.dsl)
            test_implied_vols(m, m.s0 .+ m.dsu)
            #
            m = pvm.calibrated_model(s_atm, σ_atm, T, dsl, dsu, dvl, dvu, rexl = nothing, rexu = nothing)
            c = pvm.call_option(m, m.s0)
            p = pvm.put_option(m, m.s0)
            #
            atm_vega = sqrt(T) * pdf(Normal(), 0.0)
            @test isapprox(c/atm_vega, p/atm_vega, atol=1.0e-10 )
            @test isapprox((c+p)/atm_vega, 2*σ_atm, atol=1.0e-10 )
            test_implied_vols(m, m.s0 .- m.dsl)
            test_implied_vols(m, m.s0 .+ m.dsu)
        end
    end

    @testset "Smile model" begin
        s_atm = 0.03
        σ_atm = 0.01
        r = 0.2
        dsl = [ 0.01, 0.02, 0.03, ]
        dsu = [ 0.01, 0.02, 0.03, ]
        dvl = r .* dsu
        dvu = r .* dsu
        for T in [ 1.0/12, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0 ]
            m = pvm.calibrated_model(s_atm, σ_atm, T, dsl, dsu, dvl, dvu, rexl = r, rexu = r)
            c = pvm.call_option(m, m.s0)
            p = pvm.put_option(m, m.s0)
            #
            atm_vega = sqrt(T) * pdf(Normal(), 0.0)
            @test isapprox(c/atm_vega, p/atm_vega, atol=1.0e-10 )
            @test isapprox((c+p)/atm_vega, 2*σ_atm, atol=1.0e-10 )
            test_implied_vols(m, m.s0 .- m.dsl)
            test_implied_vols(m, m.s0 .+ m.dsu)
            #
            m = pvm.calibrated_model(s_atm, σ_atm, T, dsl, dsu, dvl, dvu, rexl = nothing, rexu = nothing)
            c = pvm.call_option(m, m.s0)
            p = pvm.put_option(m, m.s0)
            #
            atm_vega = sqrt(T) * pdf(Normal(), 0.0)
            @test isapprox(c/atm_vega, p/atm_vega, atol=1.0e-10 )
            @test isapprox((c+p)/atm_vega, 2*σ_atm, atol=1.0e-10 )
            test_implied_vols(m, m.s0 .- m.dsl)
            test_implied_vols(m, m.s0 .+ m.dsu)
        end
    end

    @testset "Smirk model" begin
        s_atm = 0.03
        σ_atm = 0.01
        r = 0.2
        dsl = [ 0.01, 0.02, 0.03, ]
        dsu = [ 0.01, 0.02, 0.03, ]
        dvl = -r .* dsu
        dvu = -r .* dsu
        for T in [ 1.0/12, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0 ]
            m = pvm.calibrated_model(s_atm, σ_atm, T, dsl, dsu, dvl, dvu, rexl = -r, rexu = -r)
            c = pvm.call_option(m, m.s0)
            p = pvm.put_option(m, m.s0)
            #
            atm_vega = sqrt(T) * pdf(Normal(), 0.0)
            @test isapprox(c/atm_vega, p/atm_vega, atol=1.0e-10 )
            @test isapprox((c+p)/atm_vega, 2*σ_atm, atol=1.0e-10 )
            test_implied_vols(m, m.s0 .- m.dsl)
            test_implied_vols(m, m.s0 .+ m.dsu)
            #
            m = pvm.calibrated_model(s_atm, σ_atm, T, dsl, dsu, dvl, dvu, rexl = nothing, rexu = nothing)
            c = pvm.call_option(m, m.s0)
            p = pvm.put_option(m, m.s0)
            #
            atm_vega = sqrt(T) * pdf(Normal(), 0.0)
            @test isapprox(c/atm_vega, p/atm_vega, atol=1.0e-10 )
            @test isapprox((c+p)/atm_vega, 2*σ_atm, atol=1.0e-10 )
            test_implied_vols(m, m.s0 .- m.dsl)
            test_implied_vols(m, m.s0 .+ m.dsu)
        end
    end

    @testset "Calibrated model from slopes" begin
        s_atm = 0.03
        σ_atm = 0.01
        std_devs_lo = [ 0.5, 1.0, 2.0, 3.0, ]
        std_devs_up = [ 0.5, 1.0, 2.0, 3.0, ]
        #
        rel_strikes = [   -0.02,   -0.01 , -0.005,  0.005,    0.01,    0.02 ]
        σ_smiles =    [  0.0090,  0.0090,  0.0095,  0.0105,  0.0110,  0.0110 ]
        α = 0.1
        lmfit_kwargs = (
            autodiff = :forwarddiff,
            maxIter  = 100,
            show_trace = false,
        )
        # for T in [ 1.0/12, ]
        for T in [ 1.0/12, 0.5, 1.0, 5.0, 10.0 ]
            (m, res) = pvm.calibrated_model_from_smile(
                s_atm, σ_atm, T,
                std_devs_lo, std_devs_up,
                rel_strikes, σ_smiles;
                rexl = nothing, rexu = nothing,
                α=α, lmfit_kwargs = lmfit_kwargs,
            )
            σ_model = [
                pvm.normal_volatility(m, strike)
                for strike in m.s0 .+ rel_strikes
            ]
            fit = (σ_model .- σ_smiles) .* 1e+4  # in bp
            @test maximum(abs.(fit)) < 3.7
            @test res.converged
            # println(maximum(abs.(fit)))
            # println("dvl:")
            # display(m.dvl)
            # println("dvu:")
            # display(m.dvu)
            # println("Fit:")
            # display(fit)
            # println("Converged: " * string(res.converged))
        end


    end

end
