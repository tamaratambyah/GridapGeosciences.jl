module LaplaceBeltramiTestsSeq
using PartitionedArrays
include("../LaplaceBeltrami.jl")

with_debug() do distribute
  LaplaceBeltrami.main(distribute,1)
end

end
