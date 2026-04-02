using Test

@testset "WaveEquation" begin include("WaveEquationTests.jl") end
@testset "TransientWave" begin include("TransientWaveEquationTests.jl") end

@testset "LinearisedShallowWater" begin include("LinearisedShallowWaterTests.jl") end
@testset "ShallowWater" begin include("ShallowWaterTests.jl") end
# include("test/Geophysical/seq/runtests.jl")
