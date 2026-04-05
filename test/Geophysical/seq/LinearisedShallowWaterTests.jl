# module LinearisedShallowWaterTestsSeq
using PartitionedArrays
include("../LinearisedShallowWater.jl")

# with_debug() do distribute
#   LinearisedShallowWater.main(distribute,1)
# end

n_ref_lvls = 4

## Serial model: 2D
models = get_refined_models(n_ref_lvls)
main(models)

# end
