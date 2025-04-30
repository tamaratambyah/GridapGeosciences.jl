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

model = UnstructuredDiscreteModel(cube_grid,topo,face_labels)
num_cell_dims(model)
num_point_dims(model)

sign_flip = FESpaces.get_sign_flip(model,cell_reffe)

Table(get_cell_permutations(topo))
get_cell_permutations(topo,1)

###
RT = RaviartThomasRefFE(Float64,QUAD,1)
V = FESpace(model, RT)


###
cell_reffe = Fill(RT,num_cells(model))
trian = Triangulation(model)
conf = DivConformity() #Conformity(cell_reffe[1])
conformity = conf
### cell_conformity
ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)

ctype_lface_own_ldofs = map(reffe->get_face_own_dofs(reffe,conf),ctype_reffe)
ctype_lface_pindex_pdofs = map(reffe->get_face_own_dofs_permutations(reffe,conf),ctype_reffe)

### get_face_own_dofs_permutations
# reffe = ctype_reffe[1]
# face_own_dofs = get_face_own_dofs(reffe,conf)
# [ [collect(Int,1:length(dofs)),]  for dofs in face_own_dofs  ]
# ctype_lface_pindex_pdofs = [ [[collect(Int,1:length(dofs)),collect(Int,1:length(dofs)),]  for dofs in face_own_dofs  ]]





D = num_dims(first(ctype_reffe))
d_ctype_num_dfaces = [ map(reffe->num_faces(get_polytope(reffe),d),ctype_reffe) for d in 0:D]
cell_conformity  = CellConformity(
    cell_ctype,
    ctype_lface_own_ldofs,
    ctype_lface_pindex_pdofs,
    d_ctype_num_dfaces)


### CellFE
# ctype_reffe, cell_ctype = compress_cell_data(cell_reffe)
ctype_num_dofs = map(num_dofs,ctype_reffe)
ctype_ldof_comp = map(reffe->get_dof_to_comp(reffe),ctype_reffe)
cell_shapefuns = get_cell_shapefuns(model,cell_reffe,conformity)
cell_dof_basis = FESpaces.get_cell_dof_basis(model,cell_reffe,conformity)
cell_shapefuns_domain = ReferenceDomain()
cell_dof_basis_domain = cell_shapefuns_domain
max_order = maximum(map(get_order,ctype_reffe))

cell_fe = CellFE(
  cell_ctype,
  ctype_num_dofs,
  ctype_ldof_comp,
  cell_conformity,
  cell_shapefuns,
  cell_dof_basis,
  cell_shapefuns_domain,
  cell_dof_basis_domain,
  max_order
)

#### compute conforming cell dofs
grid_topology = get_grid_topology(model)

cell_to_ctype = cell_fe.cell_ctype
ctype_to_ldof_to_comp = cell_fe.ctype_ldof_comp
ctype_to_num_dofs = cell_fe.ctype_num_dofs

d_to_ctype_to_ldface_to_own_ldofs = cell_conformity.d_ctype_ldface_own_ldofs
ctype_to_lface_to_own_ldofs = cell_conformity.ctype_lface_own_ldofs
ctype_to_lface_to_pindex_to_pdofs = cell_conformity.ctype_lface_pindex_pdofs


D = num_cell_dims(grid_topology)
n_faces = num_faces(grid_topology)
d_to_cell_to_dfaces = [ Table(get_faces(grid_topology,D,d)) for d in 0:D]
d_to_dface_to_cells = [ Table(get_faces(grid_topology,d,D)) for d in 0:D]
d_to_offset = get_offsets(grid_topology)

####
face_to_own_dofs, ntotal, d_to_dface_to_cell, d_to_dface_to_ldface =  FESpaces._generate_face_to_own_dofs(
  n_faces,
  cell_to_ctype,
  d_to_cell_to_dfaces,
  d_to_dface_to_cells,
  d_to_offset,
  d_to_ctype_to_ldface_to_own_ldofs)


cell_to_faces = Table(get_cell_faces(grid_topology))


cell_to_lface_to_pindex = Table(get_cell_permutations(grid_topology))
# cell_to_lface_to_pindex.data .= 1
# cell_to_lface_to_pindex

#####
cell_dofs = FESpaces.CellDofsNonOriented(
  cell_to_faces,
  cell_to_lface_to_pindex,
  cell_to_ctype,
  ctype_to_lface_to_own_ldofs,
  ctype_to_num_dofs,
  face_to_own_dofs,
  ctype_to_lface_to_pindex_to_pdofs)



# get index
cell = 6
a = cell_dofs
cache = array_cache(a)

ctype = a.cell_to_ctype[cell]
n_dofs = a.ctype_to_num_dofs[ctype]
setsize!(cache,(n_dofs,))
dofs = cache.array
lface_to_own_ldofs = a.ctype_to_lface_to_own_ldofs[ctype]
p = a.cell_to_faces.ptrs[cell]-1
for (lface, own_ldofs) in enumerate(lface_to_own_ldofs)
  face = a.cell_to_faces.data[p+lface]
  pindex = a.cell_to_lface_to_pindex.data[p+lface]

  # println(ctype)
  # println(lface)
  if pindex == 2
    println("pindex 2!")
  end

  pdofs = a.ctype_to_lface_to_pindex_to_pdofs[ctype][lface][pindex]

  q = a.face_to_own_dofs.ptrs[face]-1
  for (i,ldof) in enumerate(own_ldofs)
    j = pdofs[i]
    dof = a.face_to_own_dofs.data[q+j]
    dofs[ldof] = dof
  end
end
