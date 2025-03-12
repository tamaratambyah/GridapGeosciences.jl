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
using Plots
include("initialise.jl")


### Coarset model
_model = cube_model_3D
# _model = ref_model


panel_ids = get_panel_ids(_model)
cell_coords = get_cell_coordinates(_model)
cmaps = get_cell_map(_model)
cell_node_ids = get_cell_node_ids(_model)

coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
coords_panel1_2D = lazy_map(BumpMap(), coords_panel1)


plot_coords(coords_panel1_2D,panel_ids,"panel1_2D","x","y") #### plot the lat lons in panel 1

################################################################################
##### Test Gnomonic  mapping 2D local Cartesian on panel 1 -> latlon
panel1_latlon = lazy_map(GnomonicMap(), coords_panel1_2D)
cache = array_cache(panel1_latlon)
bm1() = lazy_collect(cache,panel1_latlon)
@benchmark bm1()


################################################################################
##### Test Sigma: standard spherical mapping from latlon <-> 3D Cartesian
panel1_latlon = lazy_map(GnomonicMap(), coords_panel1_2D)
panel1_sphere = lazy_map(Sigma(),panel1_latlon)
panelp_sphere = lazy_map(InvPanelMap(), panel1_sphere, panel_ids)
cache = array_cache(panelp_sphere)
bm2() = lazy_collect(cache,panelp_sphere)
@benchmark bm2()

mapped_grid = make_grid(get_grid_topology(_model),panelp_sphere)
writevtk(mapped_grid,dir*"/CS",append=false)


panelp_latlon = lazy_map(Sigma(), panelp_sphere)

mapped_grid = make_grid(get_grid_topology(_model),panelp_latlon)
writevtk(mapped_grid,dir*"/CS_latlon",append=false)


plot_coords(panel1_latlon,panel_ids,"panel1_latlon","longitude","latitude") #### plot the lat lons in panel 1
plot_coords(panelp_latlon,panel_ids,"latlon","longitude","latitude") #### plot the lat lons in panel i=1,..,6
plot_coords(panelp_sphere,panel_ids,"sphere","","") #### plot the lat lons in panel i=1,..,6


################################################################################
##### Test central angle  mapping 2D local Cartesian on panel 1 -> central angle
central_angles = lazy_map(CentralAngleMap(), coords_panel1_2D)
cache = array_cache(central_angles)
bm1() = lazy_collect(cache,central_angles)
@benchmark bm1()


plot_coords(central_angles,panel_ids,"panel1_central_angles","alpha","beta")


################################################################################
##### Test inverse central angle  mapping central angel -> panel 1
_coords_panel1_2D = lazy_map(InverseCentralAngleMap(), central_angles)
cache = array_cache(_coords_panel1_2D)
bm1() = lazy_collect(cache,_coords_panel1_2D)
@benchmark bm1()

plot_coords(_coords_panel1_2D,panel_ids,"panel1_coords_inv","x","y")


################################################################################
