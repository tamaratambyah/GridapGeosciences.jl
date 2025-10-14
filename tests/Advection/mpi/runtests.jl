using Test

@testset "MPIAdvectionDGUpwinding" begin include("AdvectionDGUpwindingTests.jl") end
@testset "MPIAdvectionSUPG" begin include("AdvectionSUPGTests.jl") end

# include("tests/Advection/mpi/runtests.jl")
