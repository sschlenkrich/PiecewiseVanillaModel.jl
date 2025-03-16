
using PiecewiseVanillaModel
using QuadGK
using Test
pvm = PiecewiseVanillaModel

@testset "Model Integrals" begin

    @testset "Normal model integral" begin
        s0 = 0.03
        v0 = 0.01
        w0 = 0.0
        T = 2.0
        dsl = [ 0.01, 0.02 ]
        dsu = [ 0.01, 0.02, 0.03 ]
        dvl = [ 0.0, 0.0 ]
        dvu = [ 0.0, 0.0, 0.0 ]
        m = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = 0.0, rexu = 0.0)
        for (l, u) in [ (0.0, 10.0), (0.0, 3.0), (1.0, 3.0), (-2.0, 1.0) ]
            kl = pvm.risk_factor(m, l)
            ku = pvm.risk_factor(m, u)
            cl = pvm.bachelier_price(kl, m.s0, v0, T, 1.0) + kl*(1.0 - pvm.cum_dist(m, kl))
            cu = pvm.bachelier_price(ku, m.s0, v0, T, 1.0) + ku*(1.0 - pvm.cum_dist(m, ku))
            # println(cl - cu)
            #
            I1 = pvm.integral_one(u, m.s0, v0, m.w0, T, 0.0)
            I0 = pvm.integral_one(l, m.s0, v0, m.w0, T, 0.0)
            # println(I1 - I0)
            @test isapprox(I1 - I0, cl - cu, atol=1.0e-10)
            #
            I1 = pvm.integral_one(u, m.s0 + m.dsu[1], v0, m.w0 + m.dwu[1], T, 0.0)
            I0 = pvm.integral_one(l, m.s0 + m.dsu[1], v0, m.w0 + m.dwu[1], T, 0.0)
            # println(I1 - I0)
            @test isapprox(I1 - I0, cl - cu, atol=1.0e-10)
            #
            I1 = pvm.integral_one(u, m.s0 - m.dsl[2], v0, m.w0 - m.dwl[2], T, 0.0)
            I0 = pvm.integral_one(l, m.s0 - m.dsl[2], v0, m.w0 - m.dwl[2], T, 0.0)
            # println(I1 - I0)
            @test isapprox(I1 - I0, cl - cu, atol=1.0e-10)
        end
    end


    @testset "Log-normal model integral" begin
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
        for (l, u) in [ (0.0, 10.0), (0.0, 3.0), (1.0, 3.0), (-2.0, 1.0) ]
            kl = pvm.risk_factor(m, l)
            ku = pvm.risk_factor(m, u)
            cl = pvm.black_price(kl+λ, m.s0+λ, r, T, 1.0) + kl*(1.0 - pvm.cum_dist(m, kl))
            cu = pvm.black_price(ku+λ, m.s0+λ, r, T, 1.0) + ku*(1.0 - pvm.cum_dist(m, ku))
            # println(cl - cu)
            #
            I1 = pvm.integral_one(u, m.s0, v0, m.w0, T, r)
            I0 = pvm.integral_one(l, m.s0, v0, m.w0, T, r)
            # println(I1 - I0)
            @test isapprox(I1 - I0, cl - cu, atol=1.0e-10)
            #
            I1 = pvm.integral_one(u, m.s0 + m.dsu[1], m.v0 + m.dvu[1], m.w0 + m.dwu[1], T, r)
            I0 = pvm.integral_one(l, m.s0 + m.dsu[1], m.v0 + m.dvu[1], m.w0 + m.dwu[1], T, r)
            # println(I1 - I0)
            @test isapprox(I1 - I0, cl - cu, atol=1.0e-10)
            #
            I1 = pvm.integral_one(u, m.s0 - m.dsl[2], m.v0 + m.dvl[2], m.w0 - m.dwl[2], T, r)
            I0 = pvm.integral_one(l, m.s0 - m.dsl[2], m.v0 + m.dvl[2], m.w0 - m.dwl[2], T, r)
            # println(I1 - I0)
            @test isapprox(I1 - I0, cl - cu, atol=1.0e-10)
        end
    end


    @testset "Log-normal model with slope extrapolation integral" begin
        s0 = 0.03
        v0 = 0.01
        T = 2.0
        r = 0.2
        w0 = 0.5*r*T
        dsl = [ 0.01, 0.02, ]
        dsu = [ 0.01, 0.02, ]
        dvl = -r .* dsl
        dvu = r .* dsu
        m = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = 0.0, rexu = r)
        #
        ν = 0.5*r^2*T
        λ = v0/r - s0
        # test only calls
        for (l, u) in [ (0.0, 10.0), (0.0, 3.0), (1.0, 3.0), (-2.0, 1.0) ]
            kl = pvm.risk_factor(m, l)
            ku = pvm.risk_factor(m, u)
            cl = pvm.black_price(kl+λ, m.s0+λ, r, T, 1.0) + kl*(1.0 - pvm.cum_dist(m, kl))
            cu = pvm.black_price(ku+λ, m.s0+λ, r, T, 1.0) + ku*(1.0 - pvm.cum_dist(m, ku))
            # println(cl - cu)
            #
            I1 = pvm.integral_one(u, m.s0, v0, m.w0, T, r)
            I0 = pvm.integral_one(l, m.s0, v0, m.w0, T, r)
            # println(I1 - I0)
            @test isapprox(I1 - I0, cl - cu, atol=1.0e-10)
            #
            I1 = pvm.integral_one(u, m.s0 + m.dsu[1], m.v0 + m.dvu[1], m.w0 + m.dwu[1], T, r)
            I0 = pvm.integral_one(l, m.s0 + m.dsu[1], m.v0 + m.dvu[1], m.w0 + m.dwu[1], T, r)
            # println(I1 - I0)
            @test isapprox(I1 - I0, cl - cu, atol=1.0e-10)
            #
            I1 = pvm.integral_one(u, m.s0 - m.dsl[2], m.v0 + m.dvl[2], m.w0 - m.dwl[2], T, r)
            I0 = pvm.integral_one(l, m.s0 - m.dsl[2], m.v0 + m.dvl[2], m.w0 - m.dwl[2], T, r)
            # println(I1 - I0)
            @test isapprox(I1 - I0, cl - cu, atol=1.0e-10)
        end
    end

    @testset "Normal model two-integral" begin
        s0 = 0.03
        v0 = 0.01
        w0 = 0.0
        T = 2.0
        dsl = [ 0.01, 0.02 ]
        dsu = [ 0.01, 0.02, 0.03 ]
        dvl = [ 0.0, 0.0 ]
        dvu = [ 0.0, 0.0, 0.0 ]
        m = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = 0.0, rexu = 0.0)
        #
        f1 = 1.0/sqrt(2*π*m.T)
        f(w) = pvm.risk_factor(m, w)^2 * exp(-w^2/2/m.T) * f1
        #
        for (l, u) in [ (0.0, 10.0), (0.0, 3.0), (1.0, 3.0), (-2.0, 1.0) ]
            I_ref = quadgk(f, l, u)[1]
            #
            I1 = pvm.integral_two(u, m.s0, v0, m.w0, T, 0.0)
            I0 = pvm.integral_two(l, m.s0, v0, m.w0, T, 0.0)
            # println(I1 - I0)
            @test isapprox(I1 - I0, I_ref, atol=1.0e-10)
            #
            I1 = pvm.integral_two(u, m.s0 + m.dsu[1], v0, m.w0 + m.dwu[1], T, 0.0)
            I0 = pvm.integral_two(l, m.s0 + m.dsu[1], v0, m.w0 + m.dwu[1], T, 0.0)
            # println(I1 - I0)
            @test isapprox(I1 - I0, I_ref, atol=1.0e-10)
            #
            I1 = pvm.integral_two(u, m.s0 - m.dsl[2], v0, m.w0 - m.dwl[2], T, 0.0)
            I0 = pvm.integral_two(l, m.s0 - m.dsl[2], v0, m.w0 - m.dwl[2], T, 0.0)
            # println(I1 - I0)
            @test isapprox(I1 - I0, I_ref, atol=1.0e-10)
        end
    end

    @testset "Log-normal model two-integral" begin
        s0 = 0.03
        v0 = 0.01
        T = 2.0
        r = 0.2
        w0 = 0.5*r*T
        dsl = [ 0.01, 0.02, ]
        dsu = [ 0.01, 0.02, ]
        dvl = -r .* dsl
        dvu = r .* dsu
        m = pvm.model(s0, v0, w0, T, dsl, dsu, dvl, dvu, rexl = 0.0, rexu = r)
        #
        f1 = 1.0/sqrt(2*π*m.T)
        f(w) = pvm.risk_factor(m, w)^2 * exp(-w^2/2/m.T) * f1
        #
        ν = 0.5*r^2*T
        λ = v0/r - s0
        # test only calls
        for (l, u) in [ (0.0, 10.0), (0.0, 3.0), (1.0, 3.0), (-2.0, 1.0) ]
            I_ref = quadgk(f, l, u)[1]
            #
            I1 = pvm.integral_two(u, m.s0, v0, m.w0, T, r)
            I0 = pvm.integral_two(l, m.s0, v0, m.w0, T, r)
            # println(I1 - I0)
            @test isapprox(I1 - I0, I_ref, atol=1.0e-10)
            #
            I1 = pvm.integral_two(u, m.s0 + m.dsu[1], m.v0 + m.dvu[1], m.w0 + m.dwu[1], T, r)
            I0 = pvm.integral_two(l, m.s0 + m.dsu[1], m.v0 + m.dvu[1], m.w0 + m.dwu[1], T, r)
            # println(I1 - I0)
            @test isapprox(I1 - I0, I_ref, atol=1.0e-10)
            #
            I1 = pvm.integral_two(u, m.s0 - m.dsl[2], m.v0 + m.dvl[2], m.w0 - m.dwl[2], T, r)
            I0 = pvm.integral_two(l, m.s0 - m.dsl[2], m.v0 + m.dvl[2], m.w0 - m.dwl[2], T, r)
            # println(I1 - I0)
            @test isapprox(I1 - I0, I_ref, atol=1.0e-10)
        end
    end

end
