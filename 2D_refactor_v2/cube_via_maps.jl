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
panel_ids = get_panel_ids(_model)
cell_coords = get_cell_coordinates(_model)
cmaps = get_cell_map(_model)
cell_node_ids = get_cell_node_ids(_model)



################################################################################
##### Test Rotation map where panel_id is input
coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
cache = array_cache(coords_panel1)
bm1() = lazy_collect(cache,coords_panel1)
@benchmark bm1()

print_lazy_collect(coords_panel1,cell_node_ids)

# rotate panel 1 -> panel p
coords_panelp = lazy_map(InvPanelMap(), coords_panel1, panel_ids)
cache = array_cache(coords_panelp)
bm2() = lazy_collect(cache,coords_panelp)
@benchmark bm2()

# print_lazy_collect(coords_panelp,cell_node_ids)
test_coords(coords_panelp,cell_coords)



################################################################################
##### Test Rotation map where panel_id is NOT input
panel = 4
pPanelMap() = PanelRotationMap(rp1[panel])
pInvPanelMap() = PanelRotationMap(r1p[panel])

coords_panel1 = lazy_map(pPanelMap(), cell_coords[panel])
coords_panelp = lazy_map(pInvPanelMap(), coords_panel1)
cache = array_cache(coords_panelp)
bm2() = lazy_collect(cache,coords_panelp)
@benchmark bm2()

test_coords(coords_panelp,cell_coords[panel])


################################################################################
##### Test  Bump map
# take panel p 3D coords, rotate to panel 1, bump to 2D, bump to 3D, rotate to panel p
# should recover original cell_coords
coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
coords_panel1_2D = lazy_map(BumpMap(), coords_panel1)
coords_panel1_3D = lazy_map(BumpMap(), coords_panel1_2D)
coords_panelp = lazy_map(InvPanelMap(), coords_panel1_3D, panel_ids)

cache = array_cache(coords_panelp)
bm() = lazy_collect(cache,coords_panelp)
@benchmark bm()

test_coords(coords_panelp,cell_coords)

# print_lazy_collect(coords_panelp,cell_node_ids)
# print_lazy_collect(coords_panel1_2D,cell_node_ids)


# make a new grid from mapped points
mapped_grid = make_grid(get_grid_topology(_model),coords_panelp)
writevtk(mapped_grid,dir*"/cube_model_3D",append=false)

mapped_grid = make_grid(get_grid_topology(_model),coords_panel1_2D)
writevtk(mapped_grid,dir*"/cube_model_ref_2D",append=false)

###############################################################################
#### test on refined model

_model = ref_model
panel_ids = get_panel_ids(_model)
cell_coords = get_cell_coordinates(_model)
cmaps = get_cell_map(_model)
cell_node_ids = get_cell_node_ids(_model)


coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
coords_panel1_2D = lazy_map(BumpMap(), coords_panel1)
coords_panel1_3D = lazy_map(BumpMap(), coords_panel1_2D)
coords_panelp = lazy_map(InvPanelMap(), coords_panel1_3D, panel_ids)

cache = array_cache(coords_panelp)
bm() = lazy_collect(cache,coords_panelp)
@benchmark bm()

test_coords(coords_panelp,cell_coords)

mapped_grid = make_grid(get_grid_topology(_model),coords_panelp)
writevtk(mapped_grid,dir*"/ref_cube_model_3D",append=false)
