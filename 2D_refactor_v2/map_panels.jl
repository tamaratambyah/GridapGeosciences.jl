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
include("maps.jl")

dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

function make_grid(topo::UnstructuredGridTopology{Dc,Dp},nodes::AbstractArray) where {Dc,Dp}
  # Dc = num_cell_dims(topo)
  cell_reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
  cell_node_ids = get_faces(topo,Dc,0)
  cell_type = get_cell_type(topo)

  Gridap.Geometry.UnstructuredGrid(nodes,cell_node_ids,cell_reffes,cell_type,OrientationStyle(topo))
end


### 2D cube model
coarse_cube_model_2D = UnstructuredDiscreteModel(cube_surface_1_cell_per_panel_2D()...)
nodes_2D = get_node_coordinates(coarse_cube_model_2D)
cell_ids_2D = get_cell_node_ids(coarse_cube_model_2D)
# writevtk(coarse_cube_model_2D,dir*"/cube_model",append=false)


### 3D cube_model
coarse_cube_model_3D = UnstructuredDiscreteModel(cube_surface_1_cell_per_panel()...)
nodes_3D = get_node_coordinates(coarse_cube_model_3D)
cell_ids_3D = get_cell_node_ids(coarse_cube_model_3D)

@test num_cells(coarse_cube_model_2D) == num_cells(coarse_cube_model_3D)
@test cell_ids_2D == cell_ids_3D

### get maps
A,B,b = bump_matrics()
_panel1_bump = panel1BumpMap(A,B,b)

### check mapping over 2D panel 1 -> 3D panel i
ref_ids = cell_ids_2D[1]
ref_nodes_2D = nodes_2D[ref_ids]

panel1_bump = Fill( Gridap.Fields.GenericField( _panel1_bump ), length(ref_nodes_2D))
ref_nodes_3D = lazy_map(evaluate,panel1_bump,ref_nodes_2D)

_panel_map = panelMap(rotate_panel_p_to_1,rotate_panel_1_to_p)

nodes = Vector{VectorValue{3,Float64}}(undef,length(nodes_2D))

# for the coarse model, panel_id = cell_id
for panel_id in 1:6
  panel_cell_ids = cell_ids_2D[panel_id]

  panel_nodes_3D = map(x->evaluate(_panel_map,x,panel_id), ref_nodes_3D)
  nodes[panel_cell_ids] = 1*panel_nodes_3D

  println(panel_nodes_3D)
  @test panel_nodes_3D == nodes_3D[cell_ids_3D[panel_id]]
end

mapped_grid = make_grid(get_grid_topology(coarse_cube_model_2D),nodes)
writevtk(mapped_grid,dir*"/cube_model",append=false)




################################################################################
# Apply 1 level of refinement
cube_model_2D = Gridap.Adaptivity.refine(coarse_cube_model_2D )
nodes_2D = get_node_coordinates(cube_model_2D)
cell_ids_2Dh = get_cell_node_ids(cube_model_2D)
# writevtk(cube_model_2Dh,dir*"/cube_modelh",append=false)


### 3D cube_model
cube_model_3D = Gridap.Adaptivity.refine(coarse_cube_model_3D)
nodes_3Dh = get_node_coordinates(cube_model_3D)
cell_ids_3Dh = get_cell_node_ids(cube_model_3D)

@test num_cells(cube_model_2D) == num_cells(cube_model_3D)
@test cell_ids_2D == cell_ids_3D

glue = cube_model_2D.glue
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
