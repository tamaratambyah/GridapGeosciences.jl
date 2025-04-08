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

manifold_model = ref_model
ambient_model = AmbientDiscreteModel(manifold_model)
panel_ids = get_panel_ids(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)
get_background_model(Ω_parametric) == get_background_model(Ω_ambient)

pts_ambient = get_cell_points(Ω_ambient)
pts_parametric = get_cell_points(Ω_parametric)



######## mappings
R = lazy_map(x-> PanelRotationField(r1p[x]), 1:6)
_R = lazy_map(x-> PanelRotationField(rp1[x]), 1:6)
γ = GnomonicField()
γinv = InvGnomonicField()
σ = SigmaField(r)

M = lazy_map(x->    R[x]∘ σ ∘ γ, 1:6)
Minv = lazy_map(x-> γinv ∘ σ ∘ _R[x], 1:6)


######## parametric -> ambient
f(x) = x[1]
cf_parametric = CellField(f,Ω_parametric)
cell_field = get_data(cf_parametric)
cell_invmap =  lazy_map(Reindex(Minv),panel_ids)

cell_field_phys = lazy_map(Broadcasting(∘),cell_field,cell_invmap)
cf_ambient = CellData.similar_cell_field(cf_parametric,cell_field_phys,Ω_ambient,PhysicalDomain() )
cf_ambient(pts_ambient)

writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)


######## ambient -> parametric
f(x) = x[1]*x[2]*x[3]
cf_ambient  = CellField(f,Ω_ambient)
cell_field = get_data(cf_ambient)
cell_invmap =  lazy_map(Reindex(M),panel_ids)

cell_field_phys = lazy_map(Broadcasting(∘),cell_field,cell_invmap)
cf_parametric = CellData.similar_cell_field(cf_ambient,cell_field_phys,Ω_parametric,PhysicalDomain() )
cf_parametric(pts_parametric)

writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)


######### surface operators
m_parametric = Metric(cubedsphere,Ω_parametric)


surf_grad = surface_gradient(cf_parametric,m_parametric)
surf_grad(pts_parametric)
surf_div = surface_divergence(cf_parametric,m_parametric)
surf_div(pts_parametric)

# surf_lap = surface_laplacian(cf_parametric,m)
_surf_div = divergence(surf_grad) ### asks for second derivative
