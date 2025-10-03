module LaplaceBeltramiTestsMPI
using MPI, PartitionedArrays
include("../LaplaceBeltrami.jl")

with_mpi() do distribute
  LaplaceBeltrami.main(distribute,6)
end

end
