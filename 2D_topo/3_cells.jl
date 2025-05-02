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
ncells = 3

## CCAM panel ordering
data = [ 1,2,5,6, 2,3,6,7, 7,8,3,4 ] # 3 horizontal cells, 1 internal edge flipped
ptr = generate_ptr(3)
cell_node_ids = Table(data,ptr)

polytopes = fill(QUAD,ncells)
cell_type = fill(1,ncells)

nodes_2d = a.* [
  VectorValue{2, Float64}(0.0, 0.0)
 VectorValue{2, Float64}(1.0, 0.0)
 VectorValue{2, Float64}(2.0, 0.0)
 VectorValue{2, Float64}(3.0, 0.0)
 VectorValue{2, Float64}(0.0, 3.0)
 VectorValue{2, Float64}(1.0, 3.0)
 VectorValue{2, Float64}(2.0, 3.0)
 VectorValue{2, Float64}(3.0, 3.0)
]

topo = UnstructuredGridTopology(nodes_2d,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
face_labels = FaceLabeling(topo)

reffes = LagrangianRefFE(Float64,QUAD,1)
cell_reffes=[reffes]

grid = Gridap.Geometry.UnstructuredGrid(nodes_2d,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

model = Geometry.GenericDiscreteModel(grid,topo,face_labels)
