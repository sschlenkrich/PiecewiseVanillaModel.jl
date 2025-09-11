module PiecewiseVanillaModel

using Distributions
using ForwardDiff
using LinearAlgebra
using LsqFit
using QuadGK
using Roots
using SpecialFunctions

include("Model.jl")

include("Helpers.jl")
include("Integrals.jl")
include("RiskFactorFunction.jl")
include("CalibratedModel.jl")
include("BachelierImpliedVolatility.jl")
include("BlackImpliedVolatility.jl")
include("ImpliedVolatility.jl")

call_option(m::Model, strike) = call_option_analytic(m::Model, strike)
put_option(m::Model, strike) = put_option_analytic(m::Model, strike)

end
