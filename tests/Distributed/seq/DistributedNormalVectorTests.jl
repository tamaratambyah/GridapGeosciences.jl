module DistributedNormalTestsSeq
using PartitionedArrays
include("../DistributedNormalVectorTests.jl")
with_debug() do distribute
  DistributedNormalTests.main(distribute,6)
end
end # module
