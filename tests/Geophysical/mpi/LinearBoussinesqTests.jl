# module LinearBoussinesqTestsMPI
using MPI, PartitionedArrays
include("../LinearBoussinesq.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

with_mpi() do distribute
  LinearBoussinesq.main(distribute,nprocs)
end

# end
