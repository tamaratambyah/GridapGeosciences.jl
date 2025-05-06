using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays
import Gridap.TensorValues: meas

a = 1.0


manifold_grid = ManifoldGrid(cube,cube_grid_3D,panel_ids)

num_point_dims(manifold_grid)
num_cell_dims(manifold_grid)


model = ManifoldDiscreteModel(manifold_grid)
num_point_dims(model)

ref_model = Adaptivity.refine(model)
num_point_dims(ref_model)
num_cell_dims(ref_model)

writevtk(get_ambient_grid(get_grid(ref_model)),dir*"/ref_grid",append=false)

get_panel_ids(ref_model)

ref_ref_model = Adaptivity.refine(ref_model)
num_point_dims(ref_ref_model)
num_cell_dims(ref_ref_model)
get_panel_ids(ref_ref_model)
writevtk(get_ambient_grid(get_grid(ref_ref_model)),dir*"/ref_ref_grid",append=false)
