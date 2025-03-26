using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.Fields
using Gridap.TensorValues
using Gridap.Helpers
using Gridap.Geometry
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools

include("geometry/cube_surface_1_cell_per_panel.jl")
include("maps/maps.jl")
include("maps/map_matrices.jl")
include("metric/surface_metric.jl")
include("surface_metric_and_op/quadrature.jl")

a = π/4
npanels = 1
_A , _B, _b = bump_matrics(a)
A_bump = TensorValue(_A)
B_bump = TensorValue(_B)
b_bump = VectorValue(_b)
BumpMap() = Panel1BumpMap(A_bump,B_bump,b_bump)

rotate_panel_p_to_1, rotate_panel_1_to_p = panel_rotations()

"""
single face
panel/cell no | cell_node_ids
      1       |   [1 2 3 4]


  node no     | Point in 3D
      1       |   a * (1,-1,-1)
      2       |   a * (1,1,-1)
      3       |   a * (1,-1,1)
      4       |   a * (1,1,1)

"""


nodes = a.* [
  # Point(1.0, -1.0, -1.0)  # node 1
  # Point(1.0, 1.0, -1.0)   # node 2
  Point(1.0, -1.0, 1.0)   # node 3
  Point(1.0, 1.0, 1.0)    # node 4
  Point(-1.0, -1.0, 1.0)  # node 5
  Point(-1.0, 1.0, 1.0)   # node 6
  # Point(-1.0, 1.0, -1.0)  # node 7
  # Point(-1.0, -1.0, -1.0) # node 8
]

## CCAM panel ordering
data = [ 1,2,3,4 ]
ptr = generate_ptr(npanels)
cell_node_ids = Table(data,ptr)

polytopes = fill(QUAD,npanels)
cell_type = fill(1,npanels)
reffes = LagrangianRefFE(Float64,QUAD,1)
cell_reffes=[reffes]

topo = UnstructuredGridTopology(nodes,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
face_labels = FaceLabeling(topo)
topo_grid = Gridap.Geometry.UnstructuredGrid(nodes,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

coarse_face_model =  UnstructuredDiscreteModel(topo_grid,topo,face_labels)

face_model = Adaptivity.refine(coarse_face_model)

topo_cell_coords = get_cell_coordinates(get_grid(face_model))

##### map to a panel

rp = TensorValue{3,3}(rotate_panel_p_to_1[2])
topo_coords_panel1 = lazy_map(PanelRotationMap(rp), topo_cell_coords)
coords_2D = lazy_map(BumpMap(), topo_coords_panel1)


###### create cell maps

node_coords = get_nodes_from_coords(get_grid(face_model),coords_2D)
# cell_coords = lazy_map(Broadcasting(Reindex(node_coords)),get_cell_node_ids(face_model))
cell_coords = coords_2D

ctype_shapefuns = map(get_shapefuns,cell_reffes)
cell_shapefuns = expand_cell_data(ctype_shapefuns,cell_type)
default_cell_map = lazy_map(linear_combination,cell_coords,cell_shapefuns)


cm = Fields.MemoArray(default_cell_map)
evaluate(cm[1],Point(0,0))

pt = get_cell_ref_coordinates(topo_grid)
cell_Jt = lazy_map(∇,cm)
cell_Jtx = lazy_map(evaluate,cell_Jt,pt)


grid = UnstructuredGrid(node_coords,get_cell_node_ids(face_model),get_reffes(face_model),
      get_cell_type(face_model),OrientationStyle(topo_grid))

model = UnstructuredDiscreteModel(grid)
Ω = Triangulation(ReferenceFE{2},model,get_face_labeling(model))
order = 4
m = MetricInfo(metric_func,Ω)
dΩg = Measure(m,Ω,order)
out = sum( integrate(1.0,dΩg))
out*6
4*π*r^2
