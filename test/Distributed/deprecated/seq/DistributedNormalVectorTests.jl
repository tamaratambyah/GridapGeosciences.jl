module DistributedNormalTestsSeq
using PartitionedArrays
include("../DistributedNormalVectorTests.jl")
with_debug() do distribute
  DistributedNormalTests.main(distribute,1)
end
end # module
