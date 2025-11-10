# module LinearisedShallowWaterTestsMPI
using MPI, PartitionedArrays
include("../LinearisedShallowWater.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

with_mpi() do distribute
  LinearisedShallowWater.main(distribute,nprocs;octree=false,threedims=true)
end

# end
