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
include("ManifoldGrid.jl")

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

plot_coords(coords,panel_ids)
plot_coords(phys_cell_coords,panel_ids;simName="test2")

latlon = lazy_map(GnomonicMap(),phys_cell_coords)
plot_coords(latlon,panel_ids)


# sphere_maps = lazy_map(SphereAmbientCellMap() , get_panel_ids(_model), get_cell_map(_model))
# cache = array_cache(sphere_maps)
# bm1() = lazy_collect(cache,sphere_maps)
# @benchmark bm1()
