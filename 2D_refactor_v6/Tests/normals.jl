import Gridap.Geometry: FaceToCellGlue, FaceCompressedVector, push_normal

panel_model = coarse_parametric_model()
panel_ids = get_panel_ids(panel_model)
cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)


topo = get_grid_topology(panel_model)
Dc = num_cell_dims(panel_model)
bgface_to_mask = collect(Bool, .!get_isboundary_face(topo,Dc-1))


################# LEFT
trian = BoundaryTriangulation(panel_model,bgface_to_mask,1)
pts = get_cell_points(trian)

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
trian = BoundaryTriangulation(panel_model,bgface_to_mask,2)

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

_e = (n_1 + n_2)(pts)./1
map(x->vector_length.(x),_e)


###### DEBUG BELOW HERE
######################### mapping
face_panel_ids = get_panel_ids(trian)
jacobian_cf = panelwise_cellfield(forward_jacobian,trian.trian,face_panel_ids)
inv_metric_cf = CellField(analytic_inv_metric,trian)


cell_q_x = get_cell_map(cell_grid)
face_cmap_minus = lazy_map(Reindex(cell_q_x),face_to_cell_glue.face_to_cell)

# Change of domain
face_s_q = face_glue.tface_to_mface_map
fm_Jt = lazy_map(∇,face_cmap_minus)
pinv_fm_Jt = lazy_map(Operation(pinvJt),fm_Jt)
face_s_invJt = lazy_map(∘,pinv_fm_Jt,face_s_q)
face_s_n = lazy_map(Broadcasting(Operation(push_normal)),face_s_invJt,face_s_nref_minus)

cell_vectors =  Fields.MemoArray(face_s_n)
n_2_2D = GenericCellField(cell_vectors,trian,ReferenceDomain())


_n_mapped = jacobian_cf ⋅ (inv_metric_cf ⋅ n_2_2D )
ff = Operation(sqrt)(  n_2_2D  ⋅ (inv_metric_cf ⋅ n_2_2D )  )
n_mapped = _n_mapped/ff

(n_mapped - n_2)(pts)./1
