using Test

@testset "AdvectionDGUpwinding" begin include("AdvectionDGUpwindingTests.jl") end
@testset "AdvectionSUPG" begin include("AdvectionSUPGTests.jl") end

# include("tests/Advection/seq/runtests.jl")
