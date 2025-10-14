module TransientAdvectionSUPGTestsMPI
using PartitionedArrays
include("../TransientAdvectionSUPG.jl")

with_debug() do distribute
  TransientAdvectionSUPG.main(distribute,6)
end

end
