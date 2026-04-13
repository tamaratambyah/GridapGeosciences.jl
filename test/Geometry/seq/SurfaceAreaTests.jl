
include("../SurfaceArea.jl")

## Serial model: 2D
n_ref_lvls = 4
models = get_refined_models(n_ref_lvls)
SurfaceArea.main(models)
