module AdvectionSUPGTestsSeq
using PartitionedArrays
include("../AdvectionSUPG.jl")

with_debug() do distribute
  AdvectionSUPG.main(distribute,1)
end

end
