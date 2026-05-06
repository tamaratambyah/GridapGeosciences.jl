include("../TransientShallowWater.jl")

## Serial model: 2D
n_ref_lvls = 3
radius = 1.0
models = get_refined_models(n_ref_lvls,radius)
TransientShallowWaterTests.main(models[1])
