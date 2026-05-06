include("../Hierarchy.jl")

## Serial model: 2D
n_ref_lvls = 3
radius = 1.0
coarse_model = coarse_parametric_model(radius)
HierarchyTest.main(coarse_model,n_ref_lvls)
