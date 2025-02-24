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

### 2D cube model
cube_grid,topo,face_labels = cube_surface_1_cell_per_panel_2D()
cube_model_2D = UnstructuredDiscreteModel(cube_grid,topo,face_labels)
nodes_2D = get_node_coordinates(cube_model_2D)
cell_ids_2D = get_cell_node_ids(cube_model_2D)
writevtk(cube_model_2D,dir*"/cube_model",append=false)


### 3D cube_model
cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()
cube_model_3D = UnstructuredDiscreteModel(cube_grid,topo,face_labels)
nodes_3D = get_node_coordinates(cube_model_3D)
cell_ids_3D = get_cell_node_ids(cube_model_3D)

@test num_cells(cube_model_2D) == num_cells(cube_model_3D)
@test cell_ids_2D == cell_ids_3D


### check mapping over 2D panel 1 -> 3D panel i
ids_on_panel1 = cell_ids_2D[1]
nodes_on_panel1_2D = nodes_2D[ids_on_panel1]
nodes_on_panel1_3D = map(x -> bump_panel1_2D_to_3D(x), nodes_on_panel1_2D)

# for the coarse model, panel_id = cell_id
for pidx in 1:6
  R_1p = rotate_panel_1_to_p[pidx]
  nodes_on_paneli_3D = map(x-> TensorValue(R_1p)⋅x, nodes_on_panel1_3D) # 3D nodes on panel i
  @test nodes_on_paneli_3D == nodes_3D[cell_ids_3D[pidx]]
end


################################################################################
# Apply 1 level of refinement
cube_model_2Dh = Gridap.Adaptivity.refine(cube_model_2D)
nodes_2Dh = get_node_coordinates(cube_model_2Dh)
cell_ids_2Dh = get_cell_node_ids(cube_model_2Dh)
writevtk(cube_model_2Dh,dir*"/cube_modelh",append=false)


### 3D cube_model
cube_model_3Dh = Gridap.Adaptivity.refine(cube_model_3D)
nodes_3Dh = get_node_coordinates(cube_model_3Dh)
cell_ids_3Dh = get_cell_node_ids(cube_model_3Dh)

@test num_cells(cube_model_2Dh) == num_cells(cube_model_3Dh)
@test cell_ids_2Dh == cell_ids_3Dh

glue = cube_model_2Dh.glue
child_ids = Gridap.Adaptivity.n2o_reindex(glue.n2o_cell_to_child_id,glue)

n2o_faces_map = glue.n2o_faces_map
Dc = length(n2o_faces_map)-1
is_refined    = Gridap.Adaptivity.select_refined_cells(n2o_faces_map[Dc+1])
o2n_faces_map = Gridap.Adaptivity.get_o2n_faces_map(n2o_faces_map[Dc+1])
n2o_panel_map = glue.n2o_faces_map[end]

cells_on_panel1h = findall(x->x==1,n2o_panel_map)

ids_on_panel1h = (cell_ids_2Dh[cells_on_panel1h]).data
nodes_on_panel1_2Dh = nodes_2Dh[ids_on_panel1h]
nodes_on_panel1_3Dh = map(x -> bump_panel1_2D_to_3D(x), nodes_on_panel1_2Dh)

nodes_on_panel1_3Dh == nodes_3Dh[ids_on_panel1h]


T = eltype(nodes_3Dh)
nodes = Vector{T}(undef,length(nodes_2Dh))

nodes[ids_on_panel1h] = 1*nodes_on_panel1_3Dh


for pidx = 1:6

  R_1p = rotate_panel_1_to_p[pidx]

  nodes_on_panel_p_3Dh = map(x-> TensorValue(R_1p)⋅x, nodes_on_panel1_3Dh) # 3D nodes on panel i

  cells_on_panel_p = findall(x->x==pidx,n2o_panel_map)
  ids_on_panel_p = (cell_ids_2Dh[cells_on_panel_p]).data

  @test nodes_on_panel_p_3Dh == nodes_3Dh[ids_on_panel_p]

  nodes[ids_on_panel_p] = 1*nodes_on_panel_p_3Dh

end

cell_node_ids = get_cell_node_ids(cube_model_2Dh)
cell_type = get_cell_type(cube_model_2Dh)
cell_reffes = [LagrangianRefFE(Float64,QUAD,1)]
grid_f = Gridap.Geometry.UnstructuredGrid(nodes,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

topo_f = get_grid_topology(cube_model_2Dh)
fl_f = get_face_labeling(cube_model_2Dh)
_cube_model_3Dh = UnstructuredDiscreteModel(grid_f,topo_f,fl_f)

writevtk(grid_f,dir*"/cube_modelh",append=false)
