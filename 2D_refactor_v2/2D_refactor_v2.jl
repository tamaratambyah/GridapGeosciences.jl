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
include("cube_surface_1_cell_per_panel_2D.jl")
include("cube_surface_1_cell_per_panel.jl")
include("panel_rotations.jl")
include("bump_panel1.jl")

dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


cube_grid,topo,face_labels = cube_surface_1_cell_per_panel_2D()
cube_model_2D = UnstructuredDiscreteModel(cube_grid,topo,face_labels)


### 3D cube_model
cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()
cube_model_3D = UnstructuredDiscreteModel(cube_grid,topo,face_labels)



N = num_cells(cube_model_2D)
dummy_nodes = get_node_coordinates(cube_model_2D)
cell_ids = get_cell_node_ids(cube_model_2D)
T = eltype(dummy_nodes)
nodes = Vector{T}(undef,N)

ids_on_panel1 = cell_ids[1]
nodes_on_panel1_2D = dummy_nodes[ids_on_panel1]

nodes_on_panel1_3D = map(x -> bump_panel1_2D_to_3D(x), nodes_on_panel1_2D)

pidx = 3
R_1p = rotate_panel_1_to_p[pidx]
nodes_on_paneli_3D = map(x-> TensorValue(R_1p)⋅x, nodes_on_panel1_3D)





# check nodes of panel p rotate back to panel 1
for panel in 1:6
  println("Panel ", panel)
  cells_on_panel = cell_ids[panel]
  for cell in cells_on_panel
    println(cell)
  end
end

# refinement
cube_modelh = Gridap.Adaptivity.refine(cube_model)
_cell_ids = get_cell_node_ids(cube_modelh)

cells_per_panel = 4
for panel in 1:6
  println("Panel ", panel)
  idx1 = 1 + (panel-1)*cells_per_panel
  idx2 = panel*cells_per_panel
  cells_on_panel = _cell_ids[idx1:idx2]
  println(cells_on_panel)
  # for cell in cells_on_panel
  #   # println(cell)
  # end
end


get_node_coordinates(cube_modelh)
