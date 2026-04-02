using Test

@testset "MPIAdvectionDGUpwinding" begin include("AdvectionDGUpwindingTests.jl") end
@testset "MPITransientAdvectionDGUpwinding" begin include("TransientAdvectionDGUpwindingTests.jl") end

@testset "MPIAdvectionSUPG" begin include("AdvectionSUPGTests.jl") end
@testset "MPITransientAdvectionSUPG" begin include("TransientAdvectionSUPGTests.jl") end
# include("tests/Advection/mpi/runtests.jl")
