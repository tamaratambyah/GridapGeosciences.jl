module TransientAdvectionDGUpwindingTestsSeq
using PartitionedArrays
include("../TransientAdvectionDGUpwinding.jl")

with_debug() do distribute
  TransientAdvectionDGUpwinding.main(distribute,1)
end

end
