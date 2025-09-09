import Gridap.Geometry: FaceToCellGlue, FaceCompressedVector, push_normal

Geometry.get_facet_normal(trian::BoundaryTriangulation,cell_geo_map::AbstractArray) = Geometry.get_facet_normal(trian,trian.glue,cell_geo_map)

function Geometry.get_facet_normal(trian::BoundaryTriangulation,face_to_cell_glue,cell_geo_map::AbstractArray)
  bgmodel = get_background_model(trian)
  D = num_cell_dims(bgmodel)
  cell_grid = get_grid(bgmodel)
  face_glue = get_glue(trian,Val(D))

  ## if geo_map provided, map the points to ambient space
  return get_mapped_facet_normal(cell_grid,face_glue,face_to_cell_glue,cell_geo_map)

end

function get_mapped_facet_normal(
  cell_grid::Grid,
  face_glue::FaceToFaceGlue,
  face_to_cell_glue::FaceToCellGlue,
  cell_geo_map::AbstractArray
)
  println("mapped normal")

  ## Reference normal
  function f(p)
    lface_to_n = get_facet_normal(p)
    lface_to_pindex_to_perm = get_face_vertex_permutations(p,num_cell_dims(p)-1)
    nlfaces = length(lface_to_n)
    lface_pindex_to_n = [ fill(lface_to_n[lface],length(lface_to_pindex_to_perm[lface])) for lface in 1:nlfaces ]
    lface_pindex_to_n
  end
  ptops = map(get_polytope,get_reffes(cell_grid)) #fill(QUAD,num_cells(cell_grid))
  ctype_lface_pindex_to_nref = map(f, ptops)
  # ctype_lface_pindex_to_nref = map(f, get_polytopes(cell_grid))

  face_to_nref = FaceCompressedVector(ctype_lface_pindex_to_nref,face_to_cell_glue)
  face_s_nref = lazy_map(constant_field,face_to_nref)

  # Inverse of the Jacobian transpose
  _cell_q_x = get_cell_map(cell_grid)
  cell_q_x = lazy_map(∘, cell_geo_map,_cell_q_x)
  cell_q_Jt = lazy_map(∇,cell_q_x)
  cell_q_invJt = lazy_map(Operation(pinvJt),cell_q_Jt)
  face_q_invJt = lazy_map(Reindex(cell_q_invJt),face_to_cell_glue.face_to_cell)

  # Change of domain
  face_s_q = face_glue.tface_to_mface_map
  face_s_invJt = lazy_map(∘,face_q_invJt,face_s_q)
  face_s_n = lazy_map(Broadcasting(Operation(push_normal)),face_s_invJt,face_s_nref)
  Fields.MemoArray(face_s_n)
end



panel_model = coarse_parametric_model()
# panel_model = Adaptivity.refine(panel_model)
panel_ids = get_panel_ids(panel_model)
cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)

# panel_model = get_model(panel_model)
################# LEFT
trian = BoundaryTriangulation(panel_model,1)
pts = get_cell_points(trian)

cell_vectors = Geometry.get_facet_normal(trian,trian.glue,cell_geo_map)
n_2_3D = GenericCellField(cell_vectors,trian,ReferenceDomain())

n_2_3D(pts)

################################################################################
n_2_2D = get_normal_vector(trian)

face_panel_ids = get_panel_ids(trian)
glue = get_glue(trian,Val(2))

face_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), face_panel_ids)
Jt = lazy_map(∇,face_geo_map)
J = lazy_map(Operation(transpose),Jt)

cell_q_x = get_cell_map(get_grid(panel_model))
_cell_q_x = lazy_map(Reindex(cell_q_x),trian.glue.face_to_cell)
cell_pts = lazy_map(∘,_cell_q_x,glue.tface_to_mface_map)
J_face = lazy_map(∘,J,cell_pts)

J_cf = GenericCellField(Fields.MemoArray(J_face),trian,ReferenceDomain())
J_cf(pts)

inv_cf = CellField(analytic_inv_metric,trian)

_n_mapped = J_cf ⋅ (inv_cf  ⋅ n_2_2D )
ff = Operation(sqrt)(  n_2_2D   ⋅ (inv_cf⋅ n_2_2D )  )
n_mapped = _n_mapped/ff
n_mapped(pts)

sum(n_mapped(pts) .≈ n_2_3D(pts))




##################################################


bgmodel = get_background_model(trian)
D = num_cell_dims(bgmodel)
cell_grid = get_grid(bgmodel)
face_glue = get_glue(trian,Val(D))
face_to_cell_glue = trian.glue

function f(p)
  lface_to_n = get_facet_normal(p)
  lface_to_pindex_to_perm = get_face_vertex_permutations(p,num_cell_dims(p)-1)
  nlfaces = length(lface_to_n)
  lface_pindex_to_n = [ fill(lface_to_n[lface],length(lface_to_pindex_to_perm[lface])) for lface in 1:nlfaces ]
  lface_pindex_to_n
end

ptops = map(get_polytope,get_reffes(cell_grid))
ctype_lface_pindex_to_nref = map(f, ptops)
face_to_nref = FaceCompressedVector(ctype_lface_pindex_to_nref,face_to_cell_glue)
face_s_nref_plus = lazy_map(constant_field,face_to_nref)


# Inverse of the Jacobian transpose
### Map to 3D
_cell_q_x = get_cell_map(cell_grid)
cell_q_x = lazy_map(∘, cell_geo_map, _cell_q_x)
cell_q_Jt = lazy_map(∇,cell_q_x)
cell_q_invJt = lazy_map(Operation(pinvJt),cell_q_Jt)
face_q_invJt = lazy_map(Reindex(cell_q_invJt),face_to_cell_glue.face_to_cell)

# Change of domain
face_s_q = face_glue.tface_to_mface_map
face_s_invJt = lazy_map(∘,face_q_invJt,face_s_q)
face_s_n = lazy_map(Broadcasting(Operation(push_normal)),face_s_invJt,face_s_nref_plus)

cell_vectors =  Fields.MemoArray(face_s_n)
n_1 = GenericCellField(cell_vectors,trian,ReferenceDomain())
n_1(pts)

################################# RIGHT
trian = BoundaryTriangulation(panel_model,2)

bgmodel = get_background_model(trian)
D = num_cell_dims(bgmodel)
cell_grid = get_grid(bgmodel)
face_glue = get_glue(trian,Val(D))
face_to_cell_glue = trian.glue


ptops = map(get_polytope,get_reffes(cell_grid)) #fill(QUAD,num_cells(cell_grid))
ctype_lface_pindex_to_nref = map(f, ptops)

face_to_nref = FaceCompressedVector(ctype_lface_pindex_to_nref,face_to_cell_glue)
face_s_nref_minus = lazy_map(constant_field,face_to_nref)

# Inverse of the Jacobian transpose
_cell_q_x = get_cell_map(cell_grid)
cell_q_x = lazy_map(∘,cell_geo_map,_cell_q_x)
cell_q_Jt = lazy_map(∇,cell_q_x)
cell_q_invJt = lazy_map(Operation(pinvJt),cell_q_Jt)
face_q_invJt = lazy_map(Reindex(cell_q_invJt),face_to_cell_glue.face_to_cell)

# Change of domain
face_s_q = face_glue.tface_to_mface_map
face_s_invJt = lazy_map(∘,face_q_invJt,face_s_q)
face_s_n = lazy_map(Broadcasting(Operation(push_normal)),face_s_invJt,face_s_nref_minus)

cell_vectors =  Fields.MemoArray(face_s_n)
n_2 = GenericCellField(cell_vectors,trian,ReferenceDomain())
n_2(pts)
n_2_3D = n_2

_e = (n_1 + n_2)(pts)./1
map(x->vector_length.(x),_e)

face_panel_ids = get_panel_ids(trian)
panel_cfs = [n_1, n_2, n_1+n_2]
labels = ["n1","n2","n"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)

face_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), face_panel_ids)
writevtk(trian,dir*"/ambient_model_skeleton",cellfields=cellfields,append=false,geo_map=face_geo_map)


###### DEBUG BELOW HERE
######################### 2D chart normal vectors
face_panel_ids = get_panel_ids(trian)
# jacobian_cf = panelwise_cellfield(forward_jacobian,trian.trian,face_panel_ids)
# inv_metric_cf = CellField(analytic_inv_metric,trian.trian)
# pinv_jacobian_cf = panelwise_cellfield(forward_pinv_jacobian,trian.trian,face_panel_ids)

# cell_q_x = get_cell_map(cell_grid)
# cell_q_Jt = lazy_map(∇,cell_q_x)
# cell_q_invJt = lazy_map(Operation(pinvJt),cell_q_Jt)
# face_q_invJt = lazy_map(Reindex(cell_q_invJt),face_to_cell_glue.face_to_cell)

# # Change of domain
# face_s_q = face_glue.tface_to_mface_map
# face_s_invJt = lazy_map(∘,face_q_invJt,face_s_q)
# face_s_n = lazy_map(Broadcasting(Operation(push_normal)),face_s_invJt,face_s_nref_minus)

# cell_vectors =  Fields.MemoArray(face_s_n)
# n_2_2D = GenericCellField(cell_vectors,trian,ReferenceDomain())

n_2_2D = get_normal_vector(trian)


### mapping via santi formula
face_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), face_panel_ids)
Jt = lazy_map(∇,face_geo_map)
J = lazy_map(Operation(transpose),Jt)

_cell_q_x = lazy_map(Reindex(cell_q_x),face_to_cell_glue.face_to_cell)
cell_pts = lazy_map(∘,_cell_q_x,face_s_q)
J_face = lazy_map(∘,J,cell_pts)

J_cf = GenericCellField(Fields.MemoArray(J_face),trian,ReferenceDomain())
J_cf(pts)

inv_cf = CellField(analytic_inv_metric,trian)

_n_mapped = J_cf ⋅ (inv_cf  ⋅ n_2_2D )
ff = Operation(sqrt)(  n_2_2D   ⋅ (inv_cf⋅ n_2_2D )  )
n_mapped = _n_mapped/ff
n_mapped(pts)

sum(n_mapped(pts) .≈ n_2_3D(pts))

panel_ids = get_panel_ids(panel_model)
cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)

cell_vectors = Geometry.get_facet_normal(trian,trian.glue,cell_geo_map)
n_2_3D = GenericCellField(cell_vectors,trian,ReferenceDomain())

n_2_3D = get_normal_vector(trian,cell_geo_map)

n_2_3D(pts)

# function CellData.get_normal_vector(trian::Triangulation,cell_geo_map::AbstractArray)
#   println("geo maps")
#   cell_normal = Geometry.get_facet_normal(trian,cell_geo_map)
#   CellData.get_normal_vector(trian, cell_normal)
# end

# function Geometry.get_facet_normal(trian::AdaptedTriangulation,cell_geo_map=nothing)
#   println("geo maps")
#   Geometry.get_facet_normal(trian.trian,cell_geo_map)
# end
