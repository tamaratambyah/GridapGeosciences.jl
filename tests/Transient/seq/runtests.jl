using Test

@testset "TransientAdvectionSUPG" begin include("TransientAdvectionSUPGTests.jl") end
@testset "TransientAdvectionDGUpwinding" begin include("TransientAdvectionDGUpwindingTests.jl") end
@testset "TransientWave" begin include("TransientWaveEquationTests.jl") end
# include("tests/Transient/seq/runtests.jl")
