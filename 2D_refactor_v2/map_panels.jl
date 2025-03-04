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
# include("maps.jl")

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

Dc = num_cell_dims(coarse_cube_model_2D)
Dp = num_cell_dims(coarse_cube_model_2D)

### get maps
A,B,b = bump_matrics()
_panel1_bump = panel1BumpMap(A,B,b)
_panel_map = panelMap(rotate_panel_p_to_1,rotate_panel_1_to_p)


### check mapping over 2D panel 1 -> 3D panel i
n2o_panel_map = collect(1:6)
ref_cell_ids =  findall(x->x==1,n2o_panel_map)
ref_node_ids = (cell_ids_2D[ref_cell_ids]).data
ref_nodes_2D = nodes_2D[ref_node_ids]
ref_nodes_3D = map(x->evaluate(_panel1_bump,x),ref_nodes_2D)


# panel1_bump = Fill( Gridap.Fields.GenericField( _panel1_bump ), length(ref_nodes_2D))
# ref_nodes_3D = lazy_map(evaluate,panel1_bump,ref_nodes_2D)

nodes = Vector{VectorValue{Dc+1,Float64}}(undef,length(nodes_2D))


for panel_id = 1:6

  panel_cell_ids =  findall(x->x==panel_id,n2o_panel_map)
  panel_node_ids = (cell_ids_2D[panel_cell_ids]).data

  panel_nodes_3D = map(x->evaluate(_panel_map,x,panel_id), ref_nodes_3D)

  @test panel_nodes_3D == nodes_3D[panel_node_ids]

  nodes[panel_node_ids] = 1*panel_nodes_3D

end


mapped_grid = make_grid(get_grid_topology(coarse_cube_model_2D),nodes)
writevtk(mapped_grid,dir*"/cube_model",append=false)




################################################################################
# Apply 1 level of refinement
cube_model_2D = Gridap.Adaptivity.refine(coarse_cube_model_2D )
nodes_2D = get_node_coordinates(cube_model_2D)
cell_ids_2D = get_cell_node_ids(cube_model_2D)
# writevtk(cube_model_2Dh,dir*"/cube_modelh",append=false)


### 3D cube_model
cube_model_3D = Gridap.Adaptivity.refine(coarse_cube_model_3D)
nodes_3D = get_node_coordinates(cube_model_3D)
cell_ids_3D = get_cell_node_ids(cube_model_3D)

@test num_cells(cube_model_2D) == num_cells(cube_model_3D)
@test cell_ids_2D == cell_ids_3D

glue = get_adaptivity_glue(cube_model_2D)
n2o_panel_map = glue.n2o_faces_map[end]


ref_cell_ids =  findall(x->x==1,n2o_panel_map)
ref_node_ids = (cell_ids_2D[ref_cell_ids]).data
ref_nodes_2D = nodes_2D[ref_node_ids]
ref_nodes_3D = map(x->evaluate(_panel1_bump,x),ref_nodes_2D)

nodes = Vector{VectorValue{Dc+1,Float64}}(undef,length(nodes_2D))

for panel_id = 1:6

  panel_cell_ids =  findall(x->x==panel_id,n2o_panel_map)
  panel_node_ids = (cell_ids_2D[panel_cell_ids]).data

  panel_nodes_3D = map(x->evaluate(_panel_map,x,panel_id), ref_nodes_3D)

  @test panel_nodes_3D == nodes_3D[panel_node_ids]

  nodes[panel_node_ids] = 1*panel_nodes_3D

end

mapped_grid = make_grid(get_grid_topology(cube_model_2D),nodes)
writevtk(mapped_grid,dir*"/cube_model",append=false)



################################################################################
# Apply 2 level of refinement
fine_cube_model_2D = Gridap.Adaptivity.refine(cube_model_2D )
nodes_2D = get_node_coordinates(fine_cube_model_2D)
cell_ids_2D = get_cell_node_ids(fine_cube_model_2D)
# writevtk(cube_model_2Dh,dir*"/cube_modelh",append=false)


### 3D cube_model
fine_cube_model_3D = Gridap.Adaptivity.refine(cube_model_3D)
nodes_3D = get_node_coordinates(fine_cube_model_3D)
cell_ids_3D = get_cell_node_ids(fine_cube_model_3D)

@test num_cells(fine_cube_model_2D) == num_cells(fine_cube_model_3D)
@test cell_ids_2D == cell_ids_3D

fine_glue = get_adaptivity_glue(fine_cube_model_2D)

fine_n2o_panel_map = fine_glue.n2o_faces_map[end]
