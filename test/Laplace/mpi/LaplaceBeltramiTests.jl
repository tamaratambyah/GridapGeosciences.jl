# module LaplaceBeltramiTestsMPI
using MPI, PartitionedArrays
include("../LaplaceBeltrami.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
# ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

with_mpi() do distribute
  LaplaceBeltrami.main(distribute,nprocs;octree=true,threedims=true)
end

# end
