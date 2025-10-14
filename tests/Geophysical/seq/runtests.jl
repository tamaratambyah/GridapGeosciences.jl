using Test

@testset "WaveEquation" begin include("WaveEquationTests.jl") end
@testset "TransientWave" begin include("TransientWaveEquationTests.jl") end

# include("tests/Geophysical/seq/runtests.jl")
