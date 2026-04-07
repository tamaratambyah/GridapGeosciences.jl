using MPI, PartitionedArrays
include("../DistributedNormalVectorTests3D.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

with_mpi() do distribute
  DistributedNormalTests3D.main(distribute,nprocs)
end
