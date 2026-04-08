
using MPI, PartitionedArrays
include("../TransientShallowWater.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

with_mpi() do distribute
  TransientShallowWaterTests.main(distribute,nprocs)
end
