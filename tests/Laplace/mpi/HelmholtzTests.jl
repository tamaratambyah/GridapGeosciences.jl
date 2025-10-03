module HelmholtzTestsMPI
using MPI, PartitionedArrays
include("../Helmholtz.jl")

with_mpi() do distribute
  Helmholtz.main(distribute,6)
end

end
