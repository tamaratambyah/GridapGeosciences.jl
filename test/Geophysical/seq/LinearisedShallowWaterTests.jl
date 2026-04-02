module LinearisedShallowWaterTestsSeq
using PartitionedArrays
include("../LinearisedShallowWater.jl")

with_debug() do distribute
  LinearisedShallowWater.main(distribute,1)
end

end
