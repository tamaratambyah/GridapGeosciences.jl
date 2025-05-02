using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using LinearAlgebra
using FillArrays

function generate_ptr(n)
  nvertices = 4
  ptr  = Vector{Int}(undef,n+1)
  ptr[1]=1
  for i=1:n
    ptr[i+1]=ptr[i]+nvertices
  end
  ptr
end

a = 1.0
npanels = 3

## CCAM panel ordering
data = [ 1,2,5,6, 2,3,6,7, 7,4,3,8 ] # 3 horizontal cells, 1 internal edge flipped
ptr = generate_ptr(3)
cell_node_ids = Table(data,ptr)

polytopes = fill(QUAD,npanels)
cell_type = fill(1,npanels)

nodes_2d = a.* [
  Point(-1.0, -1.0)  # node 1
  Point(1.0, -1.0)   # node 2
  Point(-1.0, 1.0)   # node 3
  Point(1.0, 1.0)    # node 4
  Point(0.0, 0.0)  # node 5
  Point(0.0, 0.0)   # node 6
  Point(0.0, 0.0)   # node 6
  Point(0.0, 0.0)   # node 6
]

topo = UnstructuredGridTopology(nodes_2d,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
face_labels = FaceLabeling(topo)

reffes = LagrangianRefFE(Float64,QUAD,1)
cell_reffes=[reffes]

cube_grid = Gridap.Geometry.UnstructuredGrid(nodes_2d,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

model = Geometry.GenericDiscreteModel(cube_grid,topo,face_labels)
