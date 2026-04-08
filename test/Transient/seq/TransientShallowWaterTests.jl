include("../TransientShallowWater.jl")

## Serial model: 2D
n_ref_lvls = 3
models = get_refined_models(n_ref_lvls)
TransientShallowWaterTests.main(models[1])
