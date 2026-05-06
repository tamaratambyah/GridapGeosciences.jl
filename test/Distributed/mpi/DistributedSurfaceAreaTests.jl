using MPI, PartitionedArrays
include("../DistributedSurfaceAreaTests.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

with_mpi() do distribute
  DistributedSurfaceAreaTests.main(distribute,nprocs)
end
