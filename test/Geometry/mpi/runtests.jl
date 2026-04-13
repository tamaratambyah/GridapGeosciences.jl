using Test
using MPI

include("../../run_mpi_tests.jl")

# MPI tests -- on 1 processor only
run_mpi_tests(@__DIR__,1)
