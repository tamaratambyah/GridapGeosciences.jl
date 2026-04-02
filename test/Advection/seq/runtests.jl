using Test

# @testset "AdvectionSUPG" begin include("AdvectionSUPGTests.jl") end
# @testset "TransientAdvectionSUPG" begin include("TransientAdvectionSUPGTests.jl") end

# @testset "AdvectionDGUpwinding" begin include("AdvectionDGUpwindingTests.jl") end
@testset "TransientAdvectionDGUpwinding" begin include("TransientAdvectionDGUpwindingTests.jl") end
# include("test/Advection/seq/runtests.jl")
