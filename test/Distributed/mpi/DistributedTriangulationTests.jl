using MPI, PartitionedArrays
include("../DistributedTriangulationTests.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

with_mpi() do distribute
  DistributedTriangulationTests.main(distribute,nprocs)
end
