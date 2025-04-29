using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays


include("../src/initialise.jl")


cube_model_3D = Geometry.GenericDiscreteModel(cube_surface_1_cell_per_panel(a)...)

######### cubed sphere model
coarse_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
# model = Adaptivity.refine(coarse_model)
# ref_model = Adaptivity.refine(Adaptivity.refine(model))

# num_point_dims(coarse_model)
grid_topology = get_grid_topology(coarse_model)
num_dims(grid_topology)
num_cell_dims(grid_topology)
num_point_dims(grid_topology)
num_faces(grid_topology,0)

get_faces(grid_topology,2,0)
get_faces(grid_topology,2,1)
get_faces(grid_topology,2,2)
get_cell_faces(grid_topology)
get_offsets(grid_topology)

get_cell_permutations(grid_topology)

W = TestFESpace(coarse_model,ReferenceFE(lagrangian,Float64,1), conformity=:H1)
V = TestFESpace(coarse_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)


cell_to_faces = Table(get_cell_faces(grid_topology))
cell_to_lface_to_pindex = Table(get_cell_permutations(grid_topology))

D = num_cell_dims(grid_topology)
n_faces = num_faces(grid_topology)
d_to_cell_to_dfaces = [ Table(get_faces(grid_topology,D,d)) for d in 0:D]
d_to_dface_to_cells = [ Table(get_faces(grid_topology,d,D)) for d in 0:D]
d_to_offset = get_offsets(grid_topology)

reffes = [RaviartThomasRefFE(Float64,QUAD,1)]
cell_reffe = expand_cell_data(reffes,get_cell_type(grid_topology))

face_to_own_dofs, ntotal, d_to_dface_to_cell, d_to_dface_to_ldface =  FESpaces._generate_face_to_own_dofs(
  n_faces,
  cell_to_ctype,
  d_to_cell_to_dfaces,
  d_to_dface_to_cells,
  d_to_offset,
  d_to_ctype_to_ldface_to_own_ldofs)
