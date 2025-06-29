using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields

model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,3,0,1),(3,1),isperiodic=(false,false)))
writevtk(model,"perioidic",append=false)

cmaps = get_cell_map(model)
evaluate(cmaps[3],Point(1,1))
Ω = Triangulation(model)
dΩ = Measure(Ω,2)
u(x) = x[1]

dc = ∫(u )dΩ
get_array(dc)[3]




grid = get_grid(model)
topo = get_grid_topology(model)
labels = get_face_labeling(model)

node_coords = get_node_coordinates(grid)
cell_coords = get_cell_coordinates(grid)./1
cell_node_ids = get_cell_node_ids(grid)

vertex_coords = get_vertex_coordinates(topo)
nodes = Geometry.get_faces(topo,2,0)
edges = Geometry.get_faces(topo,2,1)
cells = Geometry.get_faces(topo,2,2)
boundaries = Geometry.get_isboundary_face(topo,1)


cell_to_entity = get_face_entity(labels,2) # For each cell, its associated entity
edge_to_entity = get_face_entity(labels,1) # For each edge, its associated entity
node_to_entity = get_face_entity(labels,0) # For each node, its associated enti

######### build by hand ########################################################
ptr = generate_ptr(3)

## topo
topo_data = [ 1,2,4,5, 2,3,5,6, 3,1,6,4 ]
_cell_vertices = Table(topo_data,ptr)
_vertex_coords = [Point(0.0,0.0),Point(1.0,0.0),Point(2.0,0.0),
                  Point(0.0,1.0),Point(1.0,1.0),Point(2.0,1.0)]
cell_type = fill(1,3)
polytopes = fill(QUAD,3)

my_topo =  UnstructuredGridTopology(_vertex_coords,_cell_vertices,cell_type,polytopes,NonOriented())
my_vertex_coords = get_vertex_coordinates(my_topo)
my_nodes = Geometry.get_faces(my_topo,2,0)
my_edges = Geometry.get_faces(my_topo,2,1)
my_cells = Geometry.get_faces(my_topo,2,2)
my_boundardies = Geometry.get_isboundary_face(my_topo,1)

my_labels = FaceLabeling(my_topo)

## grid
node_data = [ 1,2,5,6, 2,3,6,7, 3,4,7,8  ] #
_cell_node_ids = Table(node_data,ptr)
_node_coords = [Point(0.0,0.0),Point(1.0,0.0),Point(2.0,0.0), Point(3.0,0.0),
                Point(0.0,1.0),Point(1.0,1.0),Point(2.0,1.0), Point(3.0,1.0)]
cell_reffes=[LagrangianRefFE(Float64,QUAD,1)]

my_grid = UnstructuredGrid(_node_coords,_cell_node_ids,cell_reffes,cell_type,NonOriented())


my_node_coords = get_node_coordinates(my_grid)
my_cell_coords = get_cell_coordinates(my_grid)./1
my_cell_node_ids = get_cell_node_ids(my_grid)



my_model = UnstructuredDiscreteModel(my_grid,my_topo,my_labels)
writevtk(my_model,"my_perioidic",append=false)



########## coarse cube
cube_grid,cube_topo,face_labels,panel_ids = coarse_cube_surface_3D(1.0)


cell_node_ids = get_cell_node_ids(cube_grid)

nodes = get_faces(cube_topo,2,0)
edges = get_faces(cube_topo,2,1)
cells = get_faces(cube_topo,2,2)




# single_cube = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1,0,1),(1,1,1)))
# grid = get_grid(single_cube)
# topo = get_grid_topology(single_cube)

# nodes = get_node_coordinates(grid)
# vertex_coords = get_vertex_coordinates(topo)
# cell_node_ids = get_cell_node_ids(grid)
# nodes = Geometry.get_faces(topo,2,0)
# edges = Geometry.get_faces(topo,2,1)
