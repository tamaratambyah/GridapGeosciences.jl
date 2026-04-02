module ShallowWaterTestsSeq
using PartitionedArrays
include("../ShallowWater.jl")

with_debug() do distribute
  ShallowWater.main(distribute,1)
end

end
