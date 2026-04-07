using MPI, PartitionedArrays
include("../DistributedPanelIdsTests.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

with_mpi() do distribute
  DistributedPanelIdsTests.main(distribute,nprocs)
end
