using MPI, PartitionedArrays
include("../HodgeLaplacian_scalar.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

with_mpi() do distribute
  HodgeLaplacianScalarTests.main(distribute,nprocs)
end
