using MPI, PartitionedArrays
include("../HodgeLaplacian_vector.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

with_mpi() do distribute
  HodgeLaplacianVectorTests.main(distribute,nprocs)
end
