module TransientAdvectionSUPGTestsSeq
using PartitionedArrays
include("../TransientAdvectionSUPG.jl")

with_debug() do distribute
  TransientAdvectionSUPG.main(distribute,1)
end

end
