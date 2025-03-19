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


cmaps = get_cell_map(manifold_grid)
cell_coords = get_cell_coordinates(manifold_grid)
test_cell_maps(cmaps,ref_cell_coords,cell_coords)


p_cmaps = get_phys_cell_map(manifold_grid)
p_cell_coords = get_phys_cell_coordinates(manifold_grid)
test_cell_maps(p_cmaps,ref_cell_coords,p_cell_coords)

writevtk(manifold_grid.ambient_grid,dir*"/CSgrid",append=false)



using Plots
panel_ids = get_panel_ids(manifold_grid)
coords = get_cell_coordinates(manifold_grid)
phys_cell_coords = get_phys_cell_coordinates(manifold_grid)
latlon = lazy_map(GnomonicMap(),phys_cell_coords)
local_xy = lazy_map(InverseCentralAngleMap(),phys_cell_coords)
latlon_p = lazy_map(Sigma(),coords)

plot_coords(coords,panel_ids;plotTitle="sphere")
plot_coords(latlon_p,panel_ids;plotTitle="sphere_latlon")
plot_coords(phys_cell_coords,ones(Int,size(panel_ids));plotTitle="panel1_central_angles")
plot_coords(latlon,ones(Int,size(panel_ids));plotTitle="panel1_latlon")
plot_coords(local_xy,ones(Int,size(panel_ids));plotTitle="panel1_local_xy")

# sphere_maps = lazy_map(SphereAmbientCellMap() , get_panel_ids(_model), get_cell_map(_model))
# cache = array_cache(sphere_maps)
# bm1() = lazy_collect(cache,sphere_maps)
# @benchmark bm1()
