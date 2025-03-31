using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../src/initialise.jl")

manifold_model = ManifoldDiscreteModel(cube_model_3D,cube)
panel_ids = get_panel_ids(manifold_model)
order = 4
Ω = Triangulation(manifold_model)
dΩ = Measure(Ω,order)

parametric_cell_coords = get_cell_coordinates((manifold_model))
ref_cell_coords = get_cell_ref_coordinates(manifold_model)
parametric_cell_map = get_cell_map(manifold_model)
test_cell_maps(parametric_cell_map,ref_cell_coords,parametric_cell_coords)



cube_model = cube_model_3D
ref_cell_coords = get_cell_ref_coordinates(cube_model)
cmaps = get_cell_map(cube_model)
cell_coords = get_cell_coordinates(cube_model)

_coords = lazy_map(PanelMap(),cell_coords,panel_ids)
_parametric_cell_coords = lazy_map(BumpMap(),_coords)

g =  BumpField(A_bump,B_bump,b_bump)
k = map(p-> g ∘ PanelRotationField(rp1[p]), panel_ids)

parametric_cmap = lazy_map(∘,k,cmaps)

test_cell_maps(parametric_cmap,ref_cell_coords,_parametric_cell_coords)
test_cell_maps(parametric_cmap,ref_cell_coords,parametric_cell_coords)

## required for quadrature
quad = dΩ.quad
cell_Jt = lazy_map(∇,parametric_cell_map)
cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)
map(x->Gridap.TensorValues.meas(x),cell_Jtx[1])
