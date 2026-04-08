
using MPI, PartitionedArrays
include("../TransientWaveEquation.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

with_mpi() do distribute
  TransientWaveEquationTests.main(distribute,nprocs)
end
