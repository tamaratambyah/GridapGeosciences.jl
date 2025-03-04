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
include("helpers.jl")
include("cube_topo/cube_surface_1_cell_per_panel.jl")
include("maps/panel_ids_from_refinement_v2.jl")
include("maps/panel_rotations.jl")
include("maps/bump_panel1.jl")
include("maps/maps.jl")


dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

cube_model_3D = UnstructuredDiscreteModel(cube_surface_1_cell_per_panel()...)

### Coarset model
model = cube_model_3D
panel_ids = get_panel_ids(model)
cell_coords = get_cell_coordinates(model)
cmaps = get_cell_map(model)
cell_node_ids = get_cell_node_ids(model)



################################################################################
##### Test Rotation map where panel_id is input
rotate_panel_p_to_1, rotate_panel_1_to_p = panel_rotations()
rp1 = map(TensorValue,rotate_panel_p_to_1)
r1p = map(TensorValue,rotate_panel_1_to_p)

# rotate panel p -> panel 1
PanelMap() = PanelRotationMap(rp1)

coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
cache = array_cache(coords_panel1)
bm1() = lazy_collect(cache,coords_panel1)
@benchmark bm1()

print_lazy_collect(coords_panel1,cell_node_ids)

# rotate panel 1 -> panel p
InvPanelMap() = PanelRotationMap(r1p)

coords_panelp = lazy_map(InvPanelMap(), coords_panel1, panel_ids)
cache = array_cache(coords_panelp)
bm2() = lazy_collect(cache,coords_panelp)
@benchmark bm2()

# print_lazy_collect(coords_panelp,cell_node_ids)
test_coords(coords_panelp,cell_coords)



################################################################################
##### Test Rotation map where panel_id is NOT input
include("maps/maps.jl")
panel = 4
PanelMap() = PanelRotationMap(rp1[panel])
InvPanelMap() = PanelRotationMap(r1p[panel])

coords_panel1 = lazy_map(PanelMap(), cell_coords[panel])
coords_panelp = lazy_map(InvPanelMap(), coords_panel1)
cache = array_cache(coords_panelp)
bm2() = lazy_collect(cache,coords_panelp)
@benchmark bm2()

print_lazy_arr(coords_panelp)


################################################################################
##### Test  Bump map
_A , _B, _b = bump_matrics()
A_bump = TensorValue(_A)
B_bump = TensorValue(_B)
b_bump = VectorValue(_b)

BumpMap() = Panel1BumpMap(A_bump,B_bump,b_bump)

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
mapped_grid = make_grid(get_grid_topology(model),coords_panelp)
writevtk(mapped_grid,dir*"/cube_model_3D",append=false)

mapped_grid = make_grid(get_grid_topology(model),coords_panel1_2D)
writevtk(mapped_grid,dir*"/cube_model_ref_2D",append=false)

###############################################################################
#### test on refined model

# model = Gridap.Adaptivity.refine(cube_model_3D)
model = Gridap.Adaptivity.refine(model)
panel_ids = get_panel_ids(model)
cell_coords = get_cell_coordinates(model)
cmaps = get_cell_map(model)
cell_node_ids = get_cell_node_ids(model)


coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
coords_panel1_2D = lazy_map(BumpMap(), coords_panel1)
coords_panel1_3D = lazy_map(BumpMap(), coords_panel1_2D)
coords_panelp = lazy_map(InvPanelMap(), coords_panel1_3D, panel_ids)

cache = array_cache(coords_panelp)
bm() = lazy_collect(cache,coords_panelp)
@benchmark bm()

test_coords(coords_panelp,cell_coords)

mapped_grid = make_grid(get_grid_topology(model),coords_panelp)
writevtk(mapped_grid,dir*"/ref_cube_model_3D",append=false)
