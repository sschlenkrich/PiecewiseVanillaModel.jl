
using PiecewiseVanillaModel
using Test
pvm = PiecewiseVanillaModel

@testset "Vanilla options" begin

    @testset "Normal model call option" begin
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
            c_tst_1 = pvm.call_put_option_black(m1, s0 + ds, 1)
            c_tst_2 = pvm.call_put_option_black(m2, s0 + ds, 1)
            c_tst_3 = pvm.call_option_analytic(m1, s0 + ds)
            c_tst_4 = pvm.call_option_analytic(m2, s0 + ds)
            c_ref = pvm.bachelier_price(s0 + ds, s0, v0, T, 1.0)
            c_num = pvm.vanilla_option_integral(m1, s0 + ds, (s)->(s-(s0 + ds)), 1.0)
            @test isapprox(c_tst_1, c_ref, atol=1.0e-10)
            @test isapprox(c_tst_2, c_ref, atol=1.0e-10)
            @test isapprox(c_tst_3, c_ref, atol=1.0e-10)
            @test isapprox(c_tst_4, c_ref, atol=1.0e-10)
            @test isapprox(c_num,   c_ref, atol=1.0e-10)
            # println(c_tst)
            # println(c_ref)
            # println(c_num)
        end
    end

    @testset "Log-normal model call option" begin
        s0 = 0.03
        v0 = 0.01
        T = 2.0
        r = 0.2
        w0 = 0.5*r*T
        dsl = [ 0.01, 0.02, 0.04, 0.05-1.0e-6 ]
        dsu = [ 0.01, 0.02, 0.03, 0.20 ]
        dvl = -r .* dsl
        dvu = r .* dsu
        m = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu)
        #
        ν = 0.5*r^2*T
        λ = v0/r - s0
        #
        for ds in [ 0.0, 0.005, 0.01, 0.022, 0.035 ]
            c_tst_1 = pvm.call_put_option_black(m, s0 + ds, 1)
            c_tst_2 = pvm.call_option_analytic(m, s0 + ds)
            c_ref = pvm.black_price(s0 + ds + λ, s0 + λ, r, T, 1.0)
            c_num = pvm.vanilla_option_integral(m, s0 + ds, (s)->(s-(s0 + ds)), 1.0)
            @test isapprox(c_tst_1, c_ref, atol=1.0e-10)
            @test isapprox(c_tst_2, c_ref, atol=1.0e-10)
            @test isapprox(c_num, c_ref, atol=1.0e-10)
            # println(c_tst)
            # println(c_ref)
            # println(c_num)
        end
    end

    @testset "Log-normal model with slope extrapolation call option" begin
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
            c_tst_1 = pvm.call_put_option_black(m1, s0 + ds, 1)
            c_tst_2 = pvm.call_put_option_black(m2, s0 + ds, 1)
            c_tst_3 = pvm.call_option_analytic(m1, s0 + ds)
            c_tst_4 = pvm.call_option_analytic(m2, s0 + ds)
            c_ref = pvm.black_price(s0 + ds + λ, s0 + λ, r, T, 1.0)
            c_num = pvm.vanilla_option_integral(m2, s0 + ds, (s)->(s-(s0 + ds)), 1.0)
            @test isapprox(c_tst_1, c_ref, atol=1.0e-10)
            @test isapprox(c_tst_2, c_ref, atol=1.0e-10)
            @test isapprox(c_tst_3, c_ref, atol=1.0e-10)
            @test isapprox(c_tst_4, c_ref, atol=1.0e-10)
            @test isapprox(c_num, c_ref, atol=1.0e-10)
            # println(c_tst)
            # println(c_ref)
            # println(c_num)
        end
    end

    @testset "Normal model put option" begin
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
            p_tst_1 = pvm.call_put_option_black(m1, s0 - ds, -1)
            p_tst_2 = pvm.call_put_option_black(m2, s0 - ds, -1)
            p_tst_3 = pvm.put_option_analytic(m1, s0 - ds)
            p_tst_4 = pvm.put_option_analytic(m2, s0 - ds)
            p_ref = pvm.bachelier_price(s0 - ds, s0, v0, T, -1.0)
            p_num = pvm.vanilla_option_integral(m2, s0 - ds, (s)->((s0 - ds)-s), -1.0)
            @test isapprox(p_tst_1, p_ref, atol=1.0e-10)
            @test isapprox(p_tst_2, p_ref, atol=1.0e-10)
            @test isapprox(p_tst_3, p_ref, atol=1.0e-10)
            @test isapprox(p_tst_4, p_ref, atol=1.0e-10)
            @test isapprox(p_num, p_ref, atol=1.0e-10)
            # println(p_tst)
            # println(p_ref)
            # println(p_num)
        end
    end

    @testset "Log-normal model put option" begin
        s0 = 0.03
        v0 = 0.01
        T = 2.0
        r = 0.2
        w0 = 0.5*r*T
        dsl = [ 0.01, 0.02, 0.04, 0.05-1.0e-6 ]
        dsu = [ 0.01, 0.02, 0.03, 0.20 ]
        dvl = -r .* dsl
        dvu = r .* dsu
        m = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu)
        #
        ν = 0.5*r^2*T
        λ = v0/r - s0
        #
        for ds in [ 0.0, 0.005, 0.01, 0.022, 0.035 ]
            p_tst_1 = pvm.call_put_option_black(m, s0 - ds, -1)
            p_tst_2 = pvm.put_option_analytic(m, s0 - ds)
            p_ref = pvm.black_price(s0 - ds + λ, s0 + λ, r, T, -1.0)
            p_num = pvm.vanilla_option_integral(m, s0 - ds, (s)->((s0 - ds)-s), -1.0)
            @test isapprox(p_tst_1, p_ref, atol=1.0e-10)
            @test isapprox(p_tst_2, p_ref, atol=1.0e-10)
            @test isapprox(p_num, p_ref, atol=1.0e-10)
            # println(p_tst)
            # println(p_ref)
            # println(p_num)
        end
    end

    @testset "Log-normal model with slope extrapolation put option" begin
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
            p_tst_1 = pvm.call_put_option_black(m1, s0 - ds, -1)
            p_tst_2 = pvm.call_put_option_black(m2, s0 - ds, -1)
            p_tst_3 = pvm.put_option_analytic(m1, s0 - ds)
            p_tst_4 = pvm.put_option_analytic(m2, s0 - ds)
            p_ref = pvm.black_price(s0 - ds + λ, s0 + λ, r, T, -1.0)
            p_num = pvm.vanilla_option_integral(m1, s0 - ds, (s)->((s0 - ds)-s), -1.0)
            @test isapprox(p_tst_1, p_ref, atol=1.0e-10)
            @test isapprox(p_tst_2, p_ref, atol=1.0e-10)
            @test isapprox(p_tst_3, p_ref, atol=1.0e-10)
            @test isapprox(p_tst_4, p_ref, atol=1.0e-10)
            @test isapprox(p_num, p_ref, atol=1.0e-10)
            # println(p_tst)
            # println(p_ref)
            # println(p_num)
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
            c_tst_1 = pvm.call_put_option_black(m, s0 + ds, 1)
            c_tst_2 = pvm.call_option_analytic(m, s0 + ds)
            p_tst_1 = pvm.call_put_option_black(m, s0 - ds, -1)
            p_tst_2 = pvm.put_option_analytic(m, s0 - ds)
            c_num = pvm.vanilla_option_integral(m, s0 + ds, (s)->(s-(s0 + ds)), 1.0)
            p_num = pvm.vanilla_option_integral(m, s0 - ds, (s)->((s0 - ds)-s), -1.0)
            @test isapprox(c_num, c_tst_1, atol=1.0e-10)
            @test isapprox(c_num, c_tst_2, atol=1.0e-10)
            @test isapprox(p_num, p_tst_1, atol=1.0e-10)
            @test isapprox(p_num, p_tst_2, atol=1.0e-10)
            # println(c_tst)
            # println(c_num)
            # println(p_tst)
            # println(p_num)    
        end
    end


    @testset "Piece-wise model option with negative slope" begin
        s0 = 0.03
        v0 = 0.01
        T = 2.0
        r = 0.2
        w0 = 0.5*r*T
        dsl = [ 0.01, 0.02, 0.04, 0.05-1.0e-6 ]
        dsu = [ 0.01, 0.02, 0.03 ]
        dvl = -r .* dsl
        dvu = [ 0.0020, 0.0010, 0.0000 ]
        m = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu)
        for ds in [ 0.0, 0.005, 0.01, 0.022, 0.035 ]
            c_tst_1 = pvm.call_put_option_black(m, s0 + ds, 1)
            c_tst_2 = pvm.call_option_analytic(m, s0 + ds)
            p_tst_1 = pvm.call_put_option_black(m, s0 - ds, -1)
            p_tst_2 = pvm.put_option_analytic(m, s0 - ds)
            c_num = pvm.vanilla_option_integral(m, s0 + ds, (s)->(s-(s0 + ds)), 1.0)
            p_num = pvm.vanilla_option_integral(m, s0 - ds, (s)->((s0 - ds)-s), -1.0)
            @test isapprox(c_num, c_tst_1, atol=1.0e-10)
            @test isapprox(c_num, c_tst_2, atol=1.0e-10)
            @test isapprox(p_num, p_tst_1, atol=1.0e-10)
            @test isapprox(p_num, p_tst_2, atol=1.0e-10)
            # println(c_tst)
            # println(c_num)
            # println(p_tst)
            # println(p_num)
        end
    end

end
