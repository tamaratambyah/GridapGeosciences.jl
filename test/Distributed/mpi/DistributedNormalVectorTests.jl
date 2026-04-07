using MPI, PartitionedArrays
include("../DistributedNormalVectorTests.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

with_mpi() do distribute
  DistributedNormalTests.main(distribute,nprocs)
end
