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

### get the cube maps
rotate_panel_p_to_1, rotate_panel_1_to_p = panel_rotations()
rp1 = map(TensorValue,rotate_panel_p_to_1)
r1p = map(TensorValue,rotate_panel_1_to_p)

_A , _B, _b = bump_matrics()
A_bump = TensorValue(_A)
B_bump = TensorValue(_B)
b_bump = VectorValue(_b)

PanelMap() = PanelRotationMap(rp1)
InvPanelMap() = PanelRotationMap(r1p)
BumpMap() = Panel1BumpMap(A_bump,B_bump,b_bump)


### Coarset model
model = cube_model_3D
model = Gridap.Adaptivity.refine(model)


panel_ids = get_panel_ids(model)
cell_coords = get_cell_coordinates(model)
cmaps = get_cell_map(model)
cell_node_ids = get_cell_node_ids(model)

coords_panel1 = lazy_map(PanelMap(), cell_coords, panel_ids)
coords_panel1_2D = lazy_map(BumpMap(), coords_panel1)

################################################################################
##### Test Gnomonic  mapping 2D local Cartesian on panel 1 -> latlon
panel1_latlon = lazy_map(GnomonicMap(), coords_panel1_2D)
cache = array_cache(panel1_latlon)
bm1() = lazy_collect(cache,panel1_latlon)
@benchmark bm1()

# print_lazy_collect(panel1_latlon,cell_node_ids)



################################################################################
##### Test Sigma: standard spherical mapping from latlon <-> 3D Cartesian
a = 1.0
r = a/sqrt(3)
Sigma() = SigmaMap(r)

panel1_latlon = lazy_map(GnomonicMap(), coords_panel1_2D)
panel1_sphere = lazy_map(Sigma(),panel1_latlon)
panelp_sphere = lazy_map(InvPanelMap(), panel1_sphere, panel_ids)
cache = array_cache(panelp_sphere)
bm2() = lazy_collect(cache,panelp_sphere)
@benchmark bm2()

mapped_grid = make_grid(get_grid_topology(model),panelp_sphere)
writevtk(mapped_grid,dir*"/CS",append=false)


panelp_latlon = lazy_map(Sigma(), panelp_sphere)

mapped_grid = make_grid(get_grid_topology(model),panelp_latlon)
writevtk(mapped_grid,dir*"/CS_latlon",append=false)


using Plots
plot_latlons(panel1_latlon,"panel1") #### plot the lat lons in panel 1
plot_latlons(panelp_latlon,"sphere") #### plot the lat lons in panel i=1,..,6



################################################################################
##### Test central angle  mapping 2D local Cartesian on panel 1 -> central angle
central_angles = lazy_map(CentralAngleMap(), coords_panel1_2D)
cache = array_cache(central_angles)
bm1() = lazy_collect(cache,central_angles)
@benchmark bm1()

using Plots
plot_latlons(central_angles,"central_angles")
