using Test

@testset "AdvectionDGUpwinding" begin include("AdvectionDGUpwindingTests.jl") end
@testset "AdvectionSUPG" begin include("AdvectionSUPGTests.jl") end

@testset "TransientAdvectionSUPG" begin include("TransientAdvectionSUPGTests.jl") end
@testset "TransientAdvectionDGUpwinding" begin include("TransientAdvectionDGUpwindingTests.jl") end
# include("tests/Advection/seq/runtests.jl")
