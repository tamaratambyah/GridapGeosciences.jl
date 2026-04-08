using MPI, PartitionedArrays
include("../LaplaceBeltrami.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

with_mpi() do distribute
  LaplaceBeltramiTests.main(distribute,nprocs)
end
