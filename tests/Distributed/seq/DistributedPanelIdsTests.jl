module DistributedPanelIdsTestsSeq
using PartitionedArrays
include("../DistributedPanelIdsTests.jl")
with_debug() do distribute
  DistributedPanelIdsTests.main(distribute,6)
end
end # module
