using Gridap
using GridapGeosciences
using Test

include("L2ProjectionTests.jl")

n_ref_lvls = 4
## Serial model: 2D
models = get_refined_models(n_ref_lvls)
L2_projection(models)
