module DistributedPanelIdsTestsSeq
using PartitionedArrays
include("../DistributedPanelIdsTests.jl")
with_mpi() do distribute
  DistributedPanelIdsTests.main(distribute,1)
end
end # module
