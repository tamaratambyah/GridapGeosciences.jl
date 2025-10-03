module MixedHelmholtzTestsMPI
using MPI, PartitionedArrays
include("../MixedHelmholtz.jl")

with_mpi() do distribute
  MixedHelmholtz.main(distribute,6)
end

end
