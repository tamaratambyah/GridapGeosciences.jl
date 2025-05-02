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
npanels = 6

## CCAM panel ordering
data = [ 1,2,3,4, 3,4,5,6, 4,2,6,7, 6,7,5,8, 7,2,8,1, 8,1,5,3  ]
ptr = generate_ptr(npanels)
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
  Point(0.0, 0.0)  # node 7
  Point(0.0, 0.0) # node 8
]

topo = UnstructuredGridTopology(nodes_2d,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
face_labels = FaceLabeling(topo)

reffes = LagrangianRefFE(Float64,QUAD,1)
cell_reffes=[reffes]

cube_grid = Gridap.Geometry.UnstructuredGrid(nodes_2d,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

model = Geometry.GenericDiscreteModel(cube_grid,topo,face_labels)

#### permutations
reffe = RaviartThomasRefFE(Float64,QUAD,2)
conf = DivConformity()

poly = get_polytope(reffe)
face_own_dofs = get_face_own_dofs(reffe,conf)
face_vertex_perms = get_face_vertex_permutations(poly)
dranges = get_dimranges(poly)

order = get_order(get_prebasis(reffe))
freffe = LagrangianRefFE(Float64,SEGMENT,order)
f_own_dofs = get_face_own_dofs(freffe)[3]
f_dof_perms = get_face_own_dofs_permutations(freffe)[3]

_f_own_dofs = get_face_own_dofs(freffe)
_f_dof_perms = get_face_own_dofs_permutations(freffe)


out = Vector{Vector{Vector{Int}}}(undef,num_faces(poly))


for d in 1:length(dranges) # d = 1 -> vertex, d = 2 -> face, d = 3 -> cell
  dr = dranges[d]
  for j in 1:length(dr)
    face = dr[j]
    println("entity = ", face)

    trivial_perm = collect(Int,1:length(face_own_dofs[face]))
    non_trivial_perm = Int[2,1,f_own_dofs[f_dof_perms[2]]...]
    num_perms = length(face_vertex_perms[face])

    println(trivial_perm, num_perms)

    if d == 2 # face
      println("face")
      if order == 0
        out[face] = Vector{Int}[trivial_perm, trivial_perm]
      else
        out[face] = Vector{Int}[trivial_perm, non_trivial_perm]
      end
    else
      out[face] = [trivial_perm for ip in 1:num_perms]
    end

  end
end

println(out)
