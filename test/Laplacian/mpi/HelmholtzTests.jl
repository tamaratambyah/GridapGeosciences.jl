using MPI, PartitionedArrays
include("../Helmholtz.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

with_mpi() do distribute
  HelmholtzTests.main(distribute,nprocs)
end
