include("../L2ProjectionTests.jl")

## Serial model: 2D
n_ref_lvls = 4
models = get_refined_models(n_ref_lvls)
L2_projection(models)
