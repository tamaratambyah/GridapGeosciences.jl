module AdvectionSUPGTestsMPI
using PartitionedArrays
include("../AdvectionSUPG.jl")

with_debug() do distribute
  AdvectionSUPG.main(distribute,6)
end

end
