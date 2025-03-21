using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools

include("initialise.jl")


### Coarset model
_model = cube_model_3D

manifold_grid = CubeGrid(ref_ref_ref_ref_model)

ref_cell_coords = get_cell_ref_coordinates(manifold_grid)


parametric_cmaps = get_cell_map(manifold_grid)
parametric_cell_coords = get_cell_coordinates(manifold_grid)
test_cell_maps(parametric_cmaps,ref_cell_coords,parametric_cell_coords)

ambient_cmaps = get_ambient_cell_map(manifold_grid)
ambient_cell_coords = get_ambient_cell_coordinates(manifold_grid)
test_cell_maps(ambient_cmaps,ref_cell_coords,ambient_cell_coords)


# test the corner nodes of the cube have touch the sphere of radius r
# i.e. length of position vector to corner nodes == r
# corner_node_ids == 1,…,8 as all models are refined from coarse
ambient_nodes = get_ambient_node_coordinates(manifold_grid)
position_vectors = map(x-> sqrt(x[1]^2 + x[2]^2 + x[3]^2), ambient_nodes )
corner_node_ids = collect(1:8)
@test sum(position_vectors[corner_node_ids] .== r) == length(corner_node_ids)


writevtk(manifold_grid.ambient_grid,dir*"/grid",append=false)

# coords = lazy_map(CubeParametricCellMap(), get_panel_ids(_model),get_cell_map(_model))
# cache = array_cache(coords)
# bm1() = lazy_collect(cache,coords)
# @benchmark bm1()
