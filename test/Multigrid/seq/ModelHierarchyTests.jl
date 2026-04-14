include("../Hierarchy.jl")

## Serial model: 2D
n_ref_lvls = 3
coarse_model = coarse_parametric_model()
HierarchyTest.main(coarse_model,n_ref_lvls)
