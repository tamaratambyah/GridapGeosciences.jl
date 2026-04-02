# module WaveEquationTestsMPI
using MPI, PartitionedArrays
include("../WaveEquation.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

with_mpi() do distribute
  WaveEquation.main(distribute,nprocs;octree=false,threedims=true)
end

# end
