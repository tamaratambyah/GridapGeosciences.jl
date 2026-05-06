using Test

@testset "TransientShallowWater" begin include("TransientShallowWaterTests.jl") end

@testset "TransientWaveEquation" begin include("TransientWaveEquationTests.jl") end
