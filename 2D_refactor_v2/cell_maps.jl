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
# model = Gridap.Adaptivity.refine(model)
panel_ids = get_panel_ids(model)
cell_coords = get_cell_coordinates(model)
cmaps = get_cell_map(model)
cell_node_ids = get_cell_node_ids(model)
ref_cell_coords = get_cell_ref_coordinates(get_grid(model))


### get cube maps
rotate_panel_p_to_1, rotate_panel_1_to_p = panel_rotations()
rp1 = map(TensorValue,rotate_panel_p_to_1)
r1p = map(TensorValue,rotate_panel_1_to_p)


MasterPanelMaps = []
for panel = 1:6
  PanelMap() = PanelRotationMap(rp1[panel])
  InvPanelMap() = PanelRotationMap(r1p[panel])

  master = Operation(PanelMap())(InvPanelMap())
  push!(MasterPanelMaps,master)
end

cell_panel_maps = [MasterPanelMaps[i] for i in panel_ids]

cell = 4
panel_id = panel_ids[cell]

ff = Operation(cell_panel_maps[panel_id])(cmaps[cell])
@test evaluate(ff,ref_cell_coords[cell]) == cell_coords[cell]

_cell_panel_maps = lazy_map(x->x, cell_panel_maps )

gg =  lazy_map(Broadcasting(∘),_cell_panel_maps, cmaps )
evaluate(gg[1],ref_cell_coords[cell])
