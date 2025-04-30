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

get_cell_permutations(topo,2)

tables = (compute_cell_permutations(top,d) for d in 0:D)

top = topo
d = 1

D = num_cell_dims(top)
cell_to_lface_to_face = Table(get_faces(top,D,d))
data = similar(cell_to_lface_to_face.data,Int8)
ptrs = cell_to_lface_to_face.ptrs
cell_to_lface_to_pindex = Table(data,ptrs)


face_to_fvertex_to_vertex = Table(get_faces(top,d,0))
face_to_ftype = get_face_type(top,d)
reffaces = get_reffaces(Polytope{d},top)
ftype_to_pindex_to_cfvertex_to_fvertex = map(get_vertex_permutations,reffaces)
cell_to_cvertex_to_vertex = Table(get_faces(top,D,0))
cell_to_ctype = get_cell_type(top)
polytopes = get_polytopes(top)
ctype_to_lface_to_cvertices = map( (p)->get_faces(p,d,0), polytopes )

for (lface,cfvertex_to_cvertex) in enumerate(lface_to_cvertices)
  println(lface)
  println(cfvertex_to_cvertex)
end



ncells = length(cell_to_lface_to_face)
  cell = 3
    ctype = cell_to_ctype[cell]
    lface_to_cvertices = ctype_to_lface_to_cvertices[ctype]
    a = cell_to_lface_to_face.ptrs[cell]-1
    c = cell_to_cvertex_to_vertex.ptrs[cell]-1

      lface = 1
      cfvertex_to_cvertex = lface_to_cvertices[lface]

      face = cell_to_lface_to_face.data[a+lface]
      ftype = face_to_ftype[face]
      b = face_to_fvertex_to_vertex.ptrs[face]-1
      pindex_to_cfvertex_to_fvertex = ftype_to_pindex_to_cfvertex_to_fvertex[ftype]
      pindexfound = false


      pindex = 1
      cfvertex_to_fvertex = pindex_to_cfvertex_to_fvertex[pindex]

        found = true
        for (cfvertex,fvertex) in enumerate(cfvertex_to_fvertex)
          vertex1 = face_to_fvertex_to_vertex.data[b+fvertex]
          cvertex = cfvertex_to_cvertex[cfvertex]
          vertex2 = cell_to_cvertex_to_vertex.data[c+cvertex]
          if vertex1 != vertex2
            found = false
            println("breaking")
            break
          end
        end





# Geometry._compute_cell_perm_indices!(
#   cell_to_lface_to_pindex,
#   cell_to_lface_to_face,
#   cell_to_cvertex_to_vertex,
#   cell_to_ctype,
#   ctype_to_lface_to_cvertices,
#   face_to_fvertex_to_vertex,
#   face_to_ftype,
#   ftype_to_pindex_to_cfvertex_to_fvertex)

# cell_to_lface_to_pindex
