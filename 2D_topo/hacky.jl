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

nodes = a.* [
  Point(1.0, -1.0, -1.0)  # node 1
  Point(1.0, 1.0, -1.0)   # node 2
  Point(1.0, -1.0, 1.0)   # node 3
  Point(1.0, 1.0, 1.0)    # node 4
  Point(-1.0, -1.0, 1.0)  # node 5
  Point(-1.0, 1.0, 1.0)   # node 6
  Point(-1.0, 1.0, -1.0)  # node 7
  Point(-1.0, -1.0, -1.0) # node 8
]

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

Table(get_cell_permutations(topo))
get_cell_permutations(topo,1)


function FESpaces.get_face_own_dofs_permutations(reffe::ReferenceFE,conf::DivConformity)
  println("my_func")
  poly = get_polytope(reffe)
  face_own_dofs = get_face_own_dofs(reffe,conf)
  face_vertex_perms = get_face_vertex_permutations(poly)
  # _trivial_face_own_dofs_permutations(face_own_dofs)
  dranges = get_dimranges(QUAD)

  order = get_order(get_prebasis(reffe))
  freffe = LagrangianRefFE(Float64,SEGMENT,order)
  f_own_dofs = get_face_own_dofs(freffe)[3]
  f_dof_perms = get_face_own_dofs_permutations(freffe)[3]

  out = Vector{Vector{Vector{Int}}}(undef,num_faces(poly))

  for (d,r) in enumerate(dranges)
    for (iface,face) in enumerate(r)
      trivial = collect(Int,1:length(face_own_dofs[face]))

      np = length(face_vertex_perms[face])
      if d-1 != 1
        out[face] = [trivial for ip in 1:np]
      else
        if order == 0
          out[face] = Vector{Int}[trivial, trivial]
        elseif order == 1
          out[face] = Vector{Int}[trivial, Int[2,1]]
        else
          out[face] = Vector{Int}[trivial, Int[2,1,f_own_dofs[f_dof_perms[2]]...]]
        end
      end
    end
  end

  println(out)
  out
  # get_face_vertex_permutations(QUAD)
end


function FESpaces._RT_face_values(p,et,order,phi)

  # Reference facet
  @assert is_simplex(p) || is_n_cube(p) "We are assuming that all n-faces of the same n-dim are the same."
  fp = Polytope{num_dims(p)-1}(p,1)

  # geomap from ref face to polytope faces
  fgeomap = _ref_face_to_faces_geomap(p,fp)

  # Nodes are integration points (for exact integration)
  # Thus, we define the integration points in the reference
  # face polytope (fips and wips). Next, we consider the
  # n-face-wise arrays of nodes in fp (constant cell array c_fips)
  # the one of the points in the polytope after applying the geopmap
  # (fcips), and the weights for these nodes (fwips, a constant cell array)
  # Nodes (fcips)
  degree = (order)*2
  fquad = Quadrature(fp,degree)
  fips = get_coordinates(fquad)
  wips = get_weights(fquad)

  c_fips, fcips, fwips = _nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

  # Moments (fmoments)
  # The RT prebasis is expressed in terms of shape function
  fshfs = MonomialBasis(et,fp,order) ##### NEED TO CHANGE THIS TO LagrangianRefFE

  # Face moments, i.e., M(Fi)_{ab} = q_RF^a(xgp_RFi^b) w_Fi^b n_Fi ⋅ ()
  fmoments = _RT_face_moments(p, fshfs, c_fips, fcips, fwips, phi)

  return fcips, fmoments

end

get_prebasis(LagrangianRefFE(Float64,SEGMENT,1))
Gridap.Polynomials.MonomialBasis(Float64,SEGMENT,1)
get_prebasis(RT)


RT = RaviartThomasRefFE(Float64,QUAD,1)
V = FESpace(model, RT,conformity=:Hdiv)

cf = CellField(x->VectorValue(0.0,0.0),Triangulation(model))
cf(get_cell_points(Triangulation(model)))
# uh = interpolate(VectorValue(0.0,0.0),V)

trian = Triangulation(model)
pts = get_cell_points(trian)

s = get_fe_dof_basis(V)
trian = get_triangulation(s)
f = CellField(x->VectorValue(0.0,0.0),trian,DomainStyle(s))
cell_vals = s(f)


uh = interpolate(cf,V)
println(uh.cell_dof_values)
writevtk(Triangulation(model),"RT",cellfields=["u"=>uh],append=false)

poly = get_polytope(RT)
get_face_vertex_permutations(poly)


face_own_dofs = FESpaces.get_face_own_dofs(RT,DivConformity())
[ [collect(Int,1:length(dofs)),]  for dofs in face_own_dofs  ]


[[collect(Int,1:length(dofs)),collect(Int,length(dofs):-1:1),]  for dofs in face_own_dofs  ]
