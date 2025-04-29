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
ref_model = Adaptivity.refine(model)
# ref_ref_model = Adaptivity.refine(ref_model)

manifold_model = ref_model
ambient_model = AmbientDiscreteModel(manifold_model)
latlon_model = LatlonDiscreteModel(manifold_model)
panel_ids = get_panel_ids(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)
Ω_latlon = Triangulation(latlon_model)

pts_ambient = get_cell_points(Ω_ambient)
pts_parametric = get_cell_points(Ω_parametric)
pts_latlon = get_cell_points(Ω_latlon)

f_parametric(x) = x[1] #(a-x[1])*(a+x[1]) #+ (a-x[2])*(a+x[2]) ## parametric -> R
cf_parametric = CellField(f_parametric,Ω_parametric)
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)

################################################################################
## change domain: parametric -> ambient
# Mathematically, f: Ωp → R. Minv: Ωa → Ωp. Then, Minv ∘ f: (Ωa → Ωp) ∘ (Ωp → R) = Ωa → R
###############################################################################
cmap_ambient = lazy_map(x-> γinv ∘ σ ∘ PanelRotationField(rp1[x]), panel_ids)  # ambient -> parametric

# This is f ∘ Minv. Why does this work? Mathematically, this is not correct
cf_mapped = lazy_map(Broadcasting(∘),get_data(cf_parametric),cmap_ambient)

cf_ambient = CellData.similar_cell_field(cf_parametric,cf_mapped,Ω_ambient,PhysicalDomain() )
cf_ambient(pts_ambient)
writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)




# This is Minv ∘ f. Mathematically this is correct.  Why does this not work?
# _cf_mapped = lazy_map(Broadcasting(∘),cmap_ambient,get_data(cf_parametric))
# _cf_ambient = CellData.similar_cell_field(cf_parametric,_cf_mapped,Ω_ambient,PhysicalDomain() )
# _cf_ambient(pts_ambient)

################################################################################
## change domain: parametric -> ambient -> latlon
###############################################################################

# cmap_latlon = lazy_map(Reindex(Linv),panel_ids)
cmap_latlon = lazy_map(x-> γinv ∘ σ ∘ PanelRotationField(rp1[x]) ∘ σ, panel_ids) # latlon_p -> ambient -> parametric

cf_mapped = lazy_map(Broadcasting(∘),get_data(cf_parametric),cmap_latlon)

cf_latlon = CellData.similar_cell_field(cf_parametric,cf_mapped,Ω_latlon,PhysicalDomain() )
# strian = Ω_parametric
# ttrian = Ω_latlon
# D = num_cell_dims(strian)
# sglue = get_glue(strian,Val(D))
# tglue = get_glue(ttrian,Val(D))
# cf_latlon = CellData.change_domain_phys_phys(cf_parametric,ttrian,sglue,tglue)

writevtk(Ω_latlon,dir*"/latlon",cellfields=["u"=>cf_latlon],append=false)
writevtk(latlon_model,dir*"/latlon_model",append=false)

using Plots
plot_coords(get_cell_coordinates(latlon_model),panel_ids)
################################################################################
## change domain: ambient -> parametric
################################################################################
# f(x) = x[1]*x[2]*x[3] #ambient -> R
# cf_ambient  = CellField(f,Ω_ambient)

# cmap_ambient_2_parametric = lazy_map(Reindex(M),panel_ids)

# cf_mapped = lazy_map(Broadcasting(∘),get_data(cf_ambient),cmap_ambient_2_parametric)
# cf_parametric = CellData.similar_cell_field(cf_ambient,cf_mapped,Ω_parametric,PhysicalDomain() )
# cf_parametric(pts_parametric)

# writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)
# writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)
