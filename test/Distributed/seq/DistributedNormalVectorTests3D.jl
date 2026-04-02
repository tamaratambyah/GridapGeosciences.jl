module DistributedNormalTests3DSeq
using PartitionedArrays
include("../DistributedNormalVectorTests3D.jl")
with_mpi() do distribute
  DistributedNormalTests3D.main(distribute,1)
end
end # module
