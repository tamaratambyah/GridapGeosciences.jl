using Test

@testset "WaveEquation" begin include("WaveEquationTests.jl") end
@testset "LinearisedShallowWater" begin include("LinearisedShallowWaterTests.jl") end
# include("test/Geophysical/seq/runtests.jl")
