module TransientAdvectionDGUpwindingTestsMPI
using PartitionedArrays
include("../TransientAdvectionDGUpwinding.jl")

with_debug() do distribute
  TransientAdvectionDGUpwinding.main(distribute,6)
end

end
