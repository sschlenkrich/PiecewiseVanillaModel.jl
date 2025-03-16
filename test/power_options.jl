
using PiecewiseVanillaModel
using Test
pvm = PiecewiseVanillaModel

@testset "Power options" begin

    @testset "Normal model power call option" begin
        s0 = 0.03
        v0 = 0.01
        w0 = 0.0
        T = 2.0
        dsl = [ 0.01, 0.02 ]
        dsu = [ 0.01, 0.02, 0.03 ]
        dvl = [ 0.0, 0.0 ]
        dvu = [ 0.0, 0.0, 0.0 ]
        m1 = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = 0.0, rexu = 0.0)
        m2 = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = nothing, rexu = nothing)
        #
        for ds in [ 0.0, 0.005, 0.01, 0.022, 0.035 ]
            c_tst_1 = pvm.power_call_option(m1, s0 + ds)
            c_tst_2 = pvm.power_call_option(m2, s0 + ds)
            c_num = pvm.vanilla_option_integral(m1, s0 + ds, (s)->(s-(s0 + ds))^2, 1.0)
            @test isapprox(c_tst_1, c_num, atol=1.0e-10)
            @test isapprox(c_tst_2, c_num, atol=1.0e-10)
        end
    end


    @testset "Log-normal model power call option" begin
        s0 = 0.03
        v0 = 0.01
        T = 2.0
        r = 0.2
        w0 = 0.5*r*T
        dsl = [ 0.01, 0.02,]
        dsu = [ 0.01, 0.02,]
        dvl = -r .* dsl
        dvu = r .* dsu
        m1 = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = 0.0, rexu = r)
        m2 = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = nothing, rexu = nothing)
        #
        ν = 0.5*r^2*T
        λ = v0/r - s0
        #
        for ds in [ 0.0, 0.005, 0.01, 0.022, 0.035 ]
            c_tst_1 = pvm.power_call_option(m1, s0 + ds)
            c_tst_2 = pvm.power_call_option(m2, s0 + ds)
            c_num = pvm.vanilla_option_integral(m2, s0 + ds, (s)->(s-(s0 + ds))^2, 1.0)
            @test isapprox(c_tst_1, c_num, atol=1.0e-10)
            @test isapprox(c_tst_2, c_num, atol=1.0e-10)
        end
    end


    @testset "Normal model power put option" begin
        s0 = 0.03
        v0 = 0.01
        w0 = 0.0
        T = 2.0
        dsl = [ 0.01, 0.02 ]
        dsu = [ 0.01, 0.02, 0.03 ]
        dvl = [ 0.0, 0.0 ]
        dvu = [ 0.0, 0.0, 0.0 ]
        m1 = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = 0.0, rexu = 0.0)
        m2 = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = nothing, rexu = nothing)
        #
        for ds in [ 0.0, 0.005, 0.01, 0.022, 0.035 ]
            p_tst_1 = pvm.power_put_option(m1, s0 - ds)
            p_tst_2 = pvm.power_put_option(m2, s0 - ds)
            p_num = pvm.vanilla_option_integral(m2, s0 - ds, (s)->((s0 - ds)-s)^2, -1.0)
            @test isapprox(p_tst_1, p_num, atol=1.0e-10)
            @test isapprox(p_tst_2, p_num, atol=1.0e-10)
        end
    end


    @testset "Log-normal model power put option" begin
        s0 = 0.03
        v0 = 0.01
        T = 2.0
        r = 0.2
        w0 = 0.5*r*T
        dsl = [ 0.01, 0.02, ]
        dsu = [ 0.01, 0.02, ]
        dvl = -r .* dsl
        dvu = r .* dsu
        m1 = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = -r, rexu = 0.0)
        m2 = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = nothing, rexu = nothing)
        #
        ν = 0.5*r^2*T
        λ = v0/r - s0
        #
        for ds in [ 0.0, 0.005, 0.01, 0.022, 0.035 ]
            p_tst_1 = pvm.power_put_option(m1, s0 - ds)
            p_tst_2 = pvm.power_put_option(m2, s0 - ds)
            p_num = pvm.vanilla_option_integral(m1, s0 - ds, (s)->((s0 - ds)-s)^2, -1.0)
            @test isapprox(p_tst_1, p_num, atol=1.0e-10)
            @test isapprox(p_tst_2, p_num, atol=1.0e-10)
        end
    end


    @testset "Piece-wise model option" begin
        s0 = 0.03
        v0 = 0.01
        T = 2.0
        r = 0.2
        w0 = 0.5*r*T
        dsl = [ 0.01, 0.02, 0.04, 0.05-1.0e-6 ]
        dsu = [ 0.01, 0.02, 0.03 ]
        dvl = -r .* dsl
        dvu = [ 0.0020, 0.0020, 0.0030 ]
        m = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu)
        for ds in [ 0.0, 0.005, 0.01, 0.022, 0.035 ]
            c_tst = pvm.power_call_option(m, s0 + ds)
            p_tst = pvm.power_put_option(m, s0 - ds)
            c_num = pvm.vanilla_option_integral(m, s0 + ds, (s)->(s-(s0 + ds))^2, 1.0)
            p_num = pvm.vanilla_option_integral(m, s0 - ds, (s)->((s0 - ds)-s)^2, -1.0)
            @test isapprox(c_num, c_tst, atol=1.0e-10)
            @test isapprox(p_num, p_tst, atol=1.0e-10)
        end
    end


    @testset "Variance calculation" begin
        s0 = 0.03
        v0 = 0.01
        T = 2.0
        # normal model
        w0 = 0.0
        dsl = [ 0.01, ]
        dsu = [ 0.01, ]
        dvl = [ 0.0, ]
        dvu = [ 0.0, ]
        m = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = nothing, rexu = nothing)
        σ = sqrt(pvm.variance(m) / m.T)
        @test isapprox(σ, v0, atol=1.0e-10)
        # println("σ: ", σ)
        #
        # lognormal model
        r = v0 / s0
        w0 = 0.5*r*T
        dsl = [ 0.01, ]
        dsu = [ 0.01, ]
        dvl = -r .* dsl
        dvu = r .* dsu
        m = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = nothing, rexu = nothing)
        σ = sqrt(pvm.variance(m) / m.T)
        σ_ln = sqrt(s0^2 * (exp(r^2*T) - 1) / T)
        @test isapprox(σ, σ_ln, atol=1.0e-10)
        # println("σ: ", σ)
        #
        # piece-wise model
        s_atm = 0.03
        σ_atm = 0.01
        r = 0.2
        dsl = [ 0.01, 0.02, 0.04,]
        dsu = [ 0.01, 0.02, 0.03,]
        dvl = -r .* dsl
        dvu = [ 0.0020, 0.0020, 0.0030 ]
        m = pvm.calibrated_model(s_atm, σ_atm, T, dsl, dsu, dvl, dvu, rexl = nothing, rexu = 0.0)
        σ = sqrt(pvm.variance(m) / m.T)
        #
        c_num = pvm.vanilla_option_integral(m, s0, (s)->(s-s0)^2, 1.0)
        p_num = pvm.vanilla_option_integral(m, s0, (s)->(s0-s)^2, -1.0)
        σ_ref = sqrt((c_num + p_num) / m.T)
        @test isapprox(σ, σ_ref, atol=1.0e-10)
    end


    @testset "Variance versus slope" begin
        s_atm = 0.03
        σ_atm = 0.01
        T = 2.0
        dsl = [ 0.01, ]
        dsu = [ 0.01, ]
        σ_ref = [
            0.009196264690092588,
            0.009308869889107913,
            0.009440387718354950,
            0.009595268864840792,
            0.009779286016209547,
            0.010000000000000002,
            0.010267412131972839,
            0.010594883427015074,
            0.011000437542460667,
            0.011508616568579561,
            0.012153134362461359,
        ]
        σ_tst = zeros(0)
        for (k, r) in enumerate(-0.5:0.1:0.5)
            dvl = r .* dsl
            dvu = r .* dsl
            m = pvm.calibrated_model(s_atm, σ_atm, T, dsl, dsu, dvl, dvu, rexl = nothing, rexu = nothing)
            σ = sqrt(pvm.variance(m) / m.T)
            @test isapprox(σ, σ_ref[k], atol=1.0e-14)
            # println("σ: ", σ)
        end
    end


end
