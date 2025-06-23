using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields

dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

a = 1
npanels = 6

nodes_3d = a.* [
  Point(1.0, -1.0, -1.0)  # node 1
  Point(1.0, 1.0, -1.0)   # node 2
  Point(1.0, -1.0, 1.0)   # node 3
  Point(1.0, 1.0, 1.0)    # node 4
  Point(-1.0, -1.0, 1.0)  # node 5
  Point(-1.0, 1.0, 1.0)   # node 6
  Point(-1.0, 1.0, -1.0)  # node 7
  Point(-1.0, -1.0, -1.0) # node 8
]

## CCAM panel ordering
data = [ 1,2,3,4, 3,4,5,6, 2,7,4,6, 8,5,7,6, 1,8,2,7, 1,3,8,5  ] # reorient + rotated

ptr = generate_ptr(npanels)
cell_node_ids = Table(data,ptr)

polytopes = fill(QUAD,npanels)
cell_type = fill(1,npanels)
reffes = LagrangianRefFE(Float64,QUAD,1)
cell_reffes=[reffes]

topo = UnstructuredGridTopology(nodes_3d,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
_labels = FaceLabeling(topo)

nfaces = [num_faces(topo,d) for d in 0:num_cell_dims(topo)]
labels = FaceLabeling(nfaces)

points = map(x->Int32(x),collect(1:8))
edges = map(x->Int32(x),collect(9:20))
cells = map(x->Int32(x),collect(21:26))
d_to_dface_to_entity = labels.d_to_dface_to_entity
d_to_dface_to_entity[1].= points
d_to_dface_to_entity[2].= edges
d_to_dface_to_entity[3].= cells


tag_to_entities = labels.tag_to_entities

ee = map(x->[Int32(x)],collect(1:25))
eee = map(x->Int32(x),collect(1:25))
push!(ee,eee)
tag_to_name = labels.tag_to_name
labels_new = FaceLabeling(d_to_dface_to_entity,ee,tag_to_name)

cube_grid = Gridap.Geometry.UnstructuredGrid(nodes_3d,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

cube_model_3D = UnstructuredDiscreteModel(cube_grid,topo,labels_new)

Ω = Triangulation(cube_model_3D)
Γ = BoundaryTriangulation(cube_model_3D)
writevtk(cube_model_3D,dir*"/CS",append=false)



_model = CartesianDiscreteModel((0,1,0,1,0,1),(1,1,1))
ll=get_face_labeling(_model)
