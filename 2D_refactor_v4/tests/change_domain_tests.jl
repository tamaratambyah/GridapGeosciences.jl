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
ambient_model = Geometry.GenericDiscreteModel(get_ambient_grid(get_grid(manifold_model)),
get_grid_topology(manifold_model),get_face_labeling(manifold_model) )


Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)
get_background_model(Ω_parametric) == get_background_model(Ω_ambient)

pts_ambient = get_cell_points(Ω_ambient)
pts_parametric = get_cell_points(Ω_parametric)

######## ambient -> parametric
f(x) = x[1]
cf_ambient = change_domain(CellField(f,Ω_ambient),ReferenceDomain())
cf_ambient(pts_ambient)

cf_parametric = GenericCellField(CellData.get_data(cf_ambient),Ω_parametric,ReferenceDomain())
cf_parametric(pts_parametric)

writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)

######## parametric -> ambient
cf_parametric = change_domain(CellField(f,Ω_parametric),ReferenceDomain())
cf_parametric(pts_parametric)

cf_ambient = GenericCellField(CellData.get_data(cf),Ω_ambient,ReferenceDomain())
cf_ambient(pts_ambient)

writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)
