module AdvectionDGUpwindingTestsMPI
using PartitionedArrays
include("../AdvectionDGUpwinding.jl")

with_debug() do distribute
  AdvectionDGUpwinding.main(distribute,6)
end

end
