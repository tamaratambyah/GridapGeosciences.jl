module DistributedTriangulationTestsSeq
using PartitionedArrays
include("../DistributedTriangulationTests.jl")
with_debug() do distribute
  DistributedTriangulationTests.main(distribute,6)
end
end # module
