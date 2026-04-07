module DistributedSurfaceAreaTestsSeq
using PartitionedArrays
include("../DistributedSurfaceAreaTests.jl")
with_mpi() do distribute
  DistributedSurfaceAreaTests.main(distribute,1)
end
end # module
