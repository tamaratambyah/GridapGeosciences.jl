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


cmaps = get_cell_map(manifold_grid)
cell_coords = get_cell_coordinates(manifold_grid)
test_cell_maps(cmaps,ref_cell_coords,cell_coords)

p_cmaps = get_phys_cell_map(manifold_grid)
p_cell_coords = get_phys_cell_coordinates(manifold_grid)
test_cell_maps(p_cmaps,ref_cell_coords,p_cell_coords)

writevtk(manifold_grid.ambient_grid,dir*"/grid",append=false)
