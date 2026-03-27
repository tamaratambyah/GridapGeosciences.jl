module DistributedTriangulationTestsSeq
using PartitionedArrays
include("../DistributedTriangulationTests.jl")
with_mpi() do distribute
  DistributedTriangulationTests.main(distribute,1)
end
end # module
