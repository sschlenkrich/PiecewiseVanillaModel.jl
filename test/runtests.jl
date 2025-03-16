using PiecewiseVanillaModel
using Test

@testset "PiecewiseVanillaModel.jl" begin

    include("risk_factor_function.jl")
    include("integrals.jl")
    include("vanilla_options.jl")
    include("power_options.jl")
    include("calibrated_model.jl")

end
