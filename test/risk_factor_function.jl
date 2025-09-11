
using PiecewiseVanillaModel
using Test
pvm = PiecewiseVanillaModel

@testset "Model setup" begin

    @testset "Normal model setup" begin
        s0 = 0.03
        v0 = 0.01
        w0 = 0.0
        T = 2.0
        dsl = [ 0.01, 0.02 ]
        dsu = [ 0.01, 0.02, 0.03 ]
        dvl = [ 0.0, 0.0 ]
        dvu = [ 0.0, 0.0, 0.0 ]
        m = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = 0.0, rexu = 0.0)
        @test isapprox(m.dwl, [ 1.0, 2.0 ], atol=1.0e-10)
        @test isapprox(m.dwu, [ 1.0, 2.0, 3.0 ], atol=1.0e-10)
        #
        s(w) = s0 + v0*w
        w_ = [ -3.0, -2.5, -2.0, -0.5, 0.0, 0.5, 1.0, 3.0, 4.0  ]
        s_ref = [ s(w) for w in w_ ]
        s_tst = [ pvm.risk_factor(m, w) for w in w_]
        @test isapprox(s_tst, s_ref, atol=1.0e-10)
        #
        w_tst = [ pvm.brownian_factor(m, s) for s in s_tst ]
        @test isapprox(w_tst, w_, atol=1.0e-10)
    end


    @testset "Log-normal model setup" begin
        s0 = 0.03
        v0 = 0.01
        T = 2.0
        r = 0.2
        w0 = 0.5*r*T
        dsl = [ 0.01, 0.02 ]
        dsu = [ 0.01, 0.02, 0.03 ]
        dvl = -r .* dsl
        dvu = r .* dsu
        m = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu)
        #
        ν = 0.5*r^2*T
        λ = v0/r - s0
        s1(w) = s0 .+ (s0+λ).*(exp.(r.*w.-ν).-1)
        #
        @test isapprox(m.s0, s1(m.w0), atol=1.0e-10)
        @test isapprox(m.s0 .- m.dsl, s1(m.w0 .- m.dwl), atol=1.0e-10)
        @test isapprox(m.s0 .+ m.dsu, s1(m.w0 .+ m.dwu), atol=1.0e-10)
        # log-normal range
        w_ = [ m.w0-m.dwl[end], -2.0, -1.0, 0.0, w0, 0.5, 1.0, 2.0, m.w0+m.dwu[end] ]
        s_ref = [ s1(w) for w in w_ ]
        s_tst = [ pvm.risk_factor(m, w) for w in w_]
        @test isapprox(s_tst, s_ref, atol=1.0e-10)
        # lower extrapolation
        s2(w) = (s0 - dsl[end]) + (v0+dvl[end]) * (w - (w0 - m.dwl[end]))
        w_ = [ -3.0, m.w0-m.dwl[end] ]
        s_ref = [ s2(w) for w in w_ ]
        s_tst = [ pvm.risk_factor(m, w) for w in w_]
        @test isapprox(s_tst, s_ref, atol=1.0e-10)
        # upper extrapolation
        s3(w) = (s0 + dsu[end]) + (v0+dvu[end]) * (w - (w0 + m.dwu[end]))
        w_ = [ m.w0+m.dwu[end], m.w0+m.dwu[end] + 1.0 ]
        s_ref = [ s3(w) for w in w_ ]
        s_tst = [ pvm.risk_factor(m, w) for w in w_]
        @test isapprox(s_tst, s_ref, atol=1.0e-10)
        #
        w_ = [ -4.0, -3.0, -2.5, -2.0, -0.5, 0.0, 0.5, 1.0, 3.0, 4.0, 5.0 ]
        s_tst = [ pvm.risk_factor(m, w) for w in w_]
        w_tst = [ pvm.brownian_factor(m, s) for s in s_tst ]
        @test isapprox(w_tst, w_, atol=1.0e-10)
    end

    @testset "Log-normal model with slope setup" begin
        s0 = 0.03
        v0 = 0.01
        T = 2.0
        r = 0.2
        w0 = 0.5*r*T
        dsl = [ 0.01, 0.02, ]
        dsu = [ 0.01, 0.02, ]
        dvl = -r .* dsl
        dvu = r .* dsu
        m = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = -r, rexu = r)
        #
        ν = 0.5*r^2*T
        λ = v0/r - s0
        s1(w) = s0 .+ (s0+λ).*(exp.(r.*w.-ν).-1)
        #
        @test isapprox(m.s0, s1(m.w0), atol=1.0e-10)
        @test isapprox(m.s0 .- m.dsl, s1(m.w0 .- m.dwl), atol=1.0e-10)
        @test isapprox(m.s0 .+ m.dsu, s1(m.w0 .+ m.dwu), atol=1.0e-10)
        # interpolation and extrapolation range
        w_ = [ -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, w0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0 ]
        s_ref = [ s1(w) for w in w_ ]
        s_tst = [ pvm.risk_factor(m, w) for w in w_]
        # println(m.w0)
        # println(m.dwl)
        # println(m.dwu)
        @test isapprox(s_tst, s_ref, atol=1.0e-10)
        #
        s_tst = [ pvm.risk_factor(m, w) for w in w_]
        w_tst = [ pvm.brownian_factor(m, s) for s in s_tst ]
        @test isapprox(w_tst, w_, atol=1.0e-10)
    end

    @testset "Local volatility and density calculation" begin
        T = 2.0
        s_atm = 0.03
        σ_atm = 0.01
        r = 0.2
        dsl = [ 0.01, 0.02, 0.04, 0.05 ]
        dsu = [ 0.01, 0.02, 0.03 ]
        dvl = -r .* dsl
        dvu = [ 0.0020, 0.0020, 0.0030 ]
        #
        m = pvm.calibrated_model(s_atm, σ_atm, T, dsl, dsu, dvl, dvu, rexl = -r, rexu = 0.0)
        for (ds_, dv_) in zip(dsl, dvl)
            lv = pvm.local_volatility(m, s_atm - ds_)
            @test isapprox(lv, m.v0 + dv_, atol=1.1e-8)
        end
        #
        for (ds_, dv_) in zip(dsu, dvu)
            lv = pvm.local_volatility(m, s_atm + ds_)
            @test isapprox(lv, m.v0 + dv_, atol=1.1e-10)
        end
        #
        s = s_atm .+ collect(-0.06:0.005:0.04)
        #
        lv = [
            pvm.local_volatility(m, s_) for s_ in s
        ]
        lv_ref = [
            -0.001712010641763676,
            -0.000712010641763676,
             0.000287989358236323,
             0.001287984358236324,
             0.002287979358236323,
             0.003287979358236322,
             0.004287979358236323,
             0.005287979358236323,
             0.006287979358236323,
             0.007287979358236323,
             0.008287979358236323,
             0.009287979358236323,
             0.010287979358236323,
             0.011287979358236323,
             0.012287979358236324,
             0.012287979358236324,
             0.012287979358236324,
             0.012787979358236324,
             0.013287979358236324,
             0.013287979358236324,
             0.013287979358236324,
        ]
        @test isapprox(lv, lv_ref, atol=1.1e-10)
        dens = [
            pvm.implied_density(m, s_) for s_ in s
        ]
        # display(dens)
        #
        m = pvm.calibrated_model(s_atm, σ_atm, T, dsl, dsu, dvl, dvu, rexl = nothing, rexu = nothing)
        lv = [
            pvm.local_volatility(m, s_) for s_ in s
        ]
        lv_ref = [
            -0.0017217416267447288,
            -0.0007217466267447291,
             0.0002782483732552705,
             0.0012782433732552712,
             0.0022782383732552706,
             0.0032782383732552697,
             0.004278238373255271,
             0.005278238373255271,
             0.006278238373255271,
             0.007278238373255271,
             0.00827823837325527,
             0.00927823837325527,
             0.01027823837325527,
             0.01127823837325527,
             0.01227823837325527,
             0.01227823837325527,
             0.01227823837325527,
             0.012778238373255271,
             0.013278238373255272,
             0.013778238373255272,
             0.014278238373255273,
        ]
        @test isapprox(lv, lv_ref, atol=1.1e-10)
        s = s_atm .+ collect(-0.06:0.002:0.04)
        dens = [
            pvm.implied_density(m, s_) for s_ in s
        ]
        # display(dens)
    end


end
