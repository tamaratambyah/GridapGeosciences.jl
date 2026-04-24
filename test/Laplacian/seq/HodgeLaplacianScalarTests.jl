include("../HodgeLaplacian_scalar.jl")
include("../../convergence_tools.jl")

## Serial model: 2D
n_ref_lvls = 4
radius = 1.0
models = get_refined_models(n_ref_lvls,radius)
HodgeLaplacianScalarTests.main(models)
