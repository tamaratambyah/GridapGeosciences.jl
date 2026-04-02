module AdvectionDGUpwindingTestsSeq
using PartitionedArrays
include("../AdvectionDGUpwinding.jl")

with_debug() do distribute
  AdvectionDGUpwinding.main(distribute,1)
end

end
