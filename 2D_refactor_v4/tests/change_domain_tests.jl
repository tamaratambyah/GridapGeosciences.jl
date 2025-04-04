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

manifold_model = coarse_model
ambient_model = AmbientDiscreteModel(coarse_model)


Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)
get_background_model(Ω_parametric) == get_background_model(Ω_ambient)

pts_ambient = get_cell_points(Ω_ambient)
pts_parametric = get_cell_points(Ω_parametric)
f(x) = x[1]

######## parametric -> ambient
cf_parametric = CellField(f,Ω_parametric)

σ =  SigmaField(r)
σn = fill(σ,6)

cell_field = get_data(cf_parametric)
cell_invmap = lazy_map(x->σ,1:6)

cell_field_phys = lazy_map(Broadcasting(∘),cell_field,cell_invmap)
cf_ambient = CellData.similar_cell_field(cf_parametric,cell_field_phys,Ω_ambient,PhysicalDomain() )
cf_ambient(pts_ambient)

writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)




cf_ambient = CellField(f,Ω_ambient)

cell_field = get_data(cf_ambient)

cell_field_phys = lazy_map(Broadcasting(∘),cell_field,cell_invmap)
cf_parametric = CellData.similar_cell_field(cf_parametric,cell_field_phys,Ω_parametric,PhysicalDomain() )
cf_parametric(pts_parametric)
gg = gradient(cf_parametric)
gg(pts_parametric)


######## ambient -> parametric via reference domain
cf_ambient  = CellField(f,Ω_ambient)
# cf_ambient = change_domain(_cf,ReferenceDomain())
cf_ambient(pts_ambient)

cf_parametric = GenericCellField(CellData.get_data(cf_ambient),Ω_parametric,ReferenceDomain())
cf_parametric(pts_parametric)

writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)

######## parametric -> ambient via reference domain
cf_parametric = change_domain(CellField(f,Ω_parametric),ReferenceDomain())
cf_parametric(pts_parametric)

cf_ambient = GenericCellField(CellData.get_data(cf_parametric),Ω_ambient,ReferenceDomain())
cf_ambient(pts_ambient)

writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)
