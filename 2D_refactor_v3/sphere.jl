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

manifold_grid = CubedSphereGrid(ref_ref_ref_model)

ref_cell_coords = get_cell_ref_coordinates(manifold_grid)

parametric_cmaps = get_cell_map(manifold_grid)
parametric_cell_coords = get_cell_coordinates(manifold_grid)
test_cell_maps(parametric_cmaps,ref_cell_coords,parametric_cell_coords)


ambient_cmaps = get_ambient_cell_map(manifold_grid)
ambient_cell_coords = get_ambient_cell_coordinates(manifold_grid)
test_cell_maps(ambient_cmaps,ref_cell_coords,ambient_cell_coords)

writevtk(manifold_grid.ambient_grid,dir*"/CSgrid",append=false)



using Plots
panel_ids = get_panel_ids(manifold_grid)
parametric_coords = get_cell_coordinates(manifold_grid)
ambient_cell_coords = get_ambient_cell_coordinates(manifold_grid)
latlon = lazy_map(GnomonicMap(),parametric_coords)
local_xy = lazy_map(InverseCentralAngleMap(),parametric_coords)
latlon_p = lazy_map(Sigma(),ambient_cell_coords)

plot_coords(ambient_cell_coords,panel_ids;plotTitle="sphere")
plot_coords(latlon_p,panel_ids;plotTitle="sphere_latlon")
plot_coords(parametric_coords,ones(Int,size(panel_ids));plotTitle="panel1_central_angles")
plot_coords(latlon,ones(Int,size(panel_ids));plotTitle="panel1_latlon")
plot_coords(local_xy,ones(Int,size(panel_ids));plotTitle="panel1_local_xy")

# sphere_maps = lazy_map(SphereAmbientCellMap() , get_panel_ids(_model), get_cell_map(_model))
# cache = array_cache(sphere_maps)
# bm1() = lazy_collect(cache,sphere_maps)
# @benchmark bm1()
