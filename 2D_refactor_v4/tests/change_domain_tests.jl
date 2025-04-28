using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../src/initialise.jl")
coarse_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
model = Adaptivity.refine(coarse_model)
# ref_model = Adaptivity.refine(model)

manifold_model = model
ambient_model = AmbientDiscreteModel(manifold_model)
panel_ids = get_panel_ids(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)
get_background_model(Ω_parametric) == get_background_model(Ω_ambient)

pts_ambient = get_cell_points(Ω_ambient)
pts_parametric = get_cell_points(Ω_parametric)


################################################################################
## change domain: parametric -> ambient
# Mathematically, f: Ω1 → R. Minv: Ω2 → Ω1. Then, Minv ∘ f: (Ω2 → Ω1) ∘ (Ω1 → R) = Ω2 → R
################################################################################
f_parametric(x) = x[1]*x[2] # parametric -> R
cf_parametric = CellField(f_parametric,Ω_parametric)

cmap_parametric_2_ambient =  lazy_map(Reindex(Minv),panel_ids)

# This is f ∘ Minv. Why does this work? Mathematically, this is not correct
cf_mapped = lazy_map(Broadcasting(∘),get_data(cf_parametric),cmap_parametric_2_ambient)

cf_ambient = CellData.similar_cell_field(cf_parametric,cf_mapped,Ω_ambient,PhysicalDomain() )
cf_ambient(pts_ambient)
writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)

# This is Minv ∘ f. Mathematically this is correct.  Why does this not work?
_cf_mapped = lazy_map(Broadcasting(∘),cmap_parametric_2_ambient,get_data(cf_parametric))
_cf_ambient = CellData.similar_cell_field(cf_parametric,_cf_mapped,Ω_ambient,PhysicalDomain() )
_cf_ambient(pts_ambient)


################################################################################
## change domain: ambient -> parametric
################################################################################
f(x) = x[1]*x[2]*x[3] #ambient -> R
cf_ambient  = CellField(f,Ω_ambient)

cmap_ambient_2_parametric = lazy_map(Reindex(M),panel_ids)

cf_mapped = lazy_map(Broadcasting(∘),get_data(cf_ambient),cmap_ambient_2_parametric)
cf_parametric = CellData.similar_cell_field(cf_ambient,cf_mapped,Ω_parametric,PhysicalDomain() )
cf_parametric(pts_parametric)

writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)


# ######### surface operators
# m_parametric = Metric(cubedsphere,Ω_parametric)


# surf_grad = surface_gradient(cf_parametric,m_parametric)
# surf_grad(pts_parametric)
# surf_div = surface_divergence(cf_parametric,m_parametric)
# surf_div(pts_parametric)

# # surf_lap = surface_laplacian(cf_parametric,m)
# _surf_div = divergence(surf_grad) ### asks for second derivative
