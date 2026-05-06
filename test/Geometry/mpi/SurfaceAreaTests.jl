using MPI, PartitionedArrays
include("../SurfaceArea.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

with_mpi() do distribute
  SurfaceArea.main(distribute,nprocs)
end
