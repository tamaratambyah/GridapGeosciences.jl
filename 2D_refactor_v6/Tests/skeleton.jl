
################################################################################




panel_model = coarse_parametric_model()
panel_ids = get_panel_ids(panel_model)
cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)

topo = get_grid_topology(panel_model)
Dc = num_cell_dims(panel_model)
bgface_to_mask = collect(Bool, .!get_isboundary_face(topo,Dc-1))

trian = BoundaryTriangulation(panel_model,bgface_to_mask,1)

bgmodel = get_background_model(trian)
D = num_cell_dims(bgmodel)
cell_grid = get_grid(bgmodel)
face_glue = get_glue(trian,Val(D))
face_to_cell_glue = trian.glue

import Gridap.Geometry: FaceToCellGlue, FaceCompressedVector, push_normal

  function f(p)
    lface_to_n = get_facet_normal(p)
    lface_to_pindex_to_perm = get_face_vertex_permutations(p,num_cell_dims(p)-1)
    nlfaces = length(lface_to_n)
    lface_pindex_to_n = [ fill(lface_to_n[lface],length(lface_to_pindex_to_perm[lface])) for lface in 1:nlfaces ]
    lface_pindex_to_n
  end

  # ctype_lface_pindex_to_nref = map(f, get_polytopes(cell_grid))
  ptops = map(get_polytope,get_reffes(cell_grid))
  ctype_lface_pindex_to_nref = map(f, ptops)

  face_to_nref = FaceCompressedVector(ctype_lface_pindex_to_nref,face_to_cell_glue)
  face_s_nref = lazy_map(constant_field,face_to_nref)
  face_s_nref_plus = face_s_nref

  # Inverse of the Jacobian transpose
  # panel_maps = get_cell_map(cell_grid)
  # cube2panel_map = lazy_map(p-> grid_maps[p], panel_ids)
  # cell_q_x = lazy_map(∘,cube2panel_map,panel_maps)

  # evaluate(cell_q_x[6],Point(1,1))

  _cell_q_x = get_cell_map(cell_grid)
  cell_q_x = lazy_map(∘, cell_geo_map, _cell_q_x)

  cell_q_Jt = lazy_map(∇,cell_q_x)
  cell_q_invJt = lazy_map(Operation(pinvJt),cell_q_Jt)
  face_q_invJt = lazy_map(Reindex(cell_q_invJt),face_to_cell_glue.face_to_cell)

  # Change of domain
  face_s_q = face_glue.tface_to_mface_map
  face_s_invJt = lazy_map(∘,face_q_invJt,face_s_q)
  face_s_n = lazy_map(Broadcasting(Operation(push_normal)),face_s_invJt,face_s_nref)
  cell_vectors =  Fields.MemoArray(face_s_n)

n_1 = GenericCellField(cell_vectors,trian,ReferenceDomain())

# lazy_map(evaluate,face_s_q,Fill([VectorValue(0.,),VectorValue(1.,)],length(face_s_nref)))
face_cmap_plus = lazy_map(Reindex(cell_q_x),face_to_cell_glue.face_to_cell)
face_face_map_plus = lazy_map(∘,face_cmap_plus,face_s_q)
c = lazy_map(evaluate,face_face_map_plus,Fill([VectorValue(0.,),VectorValue(1.,)],length(face_s_nref)))

#################################

##### right
trian = BoundaryTriangulation(panel_model,bgface_to_mask,2)

bgmodel = get_background_model(trian)
D = num_cell_dims(bgmodel)
cell_grid = get_grid(bgmodel)
face_glue = get_glue(trian,Val(D))
face_to_cell_glue = trian.glue


  ptops = map(get_polytope,get_reffes(cell_grid)) #fill(QUAD,num_cells(cell_grid))
  ctype_lface_pindex_to_nref = map(f, ptops)

  face_to_nref = FaceCompressedVector(ctype_lface_pindex_to_nref,face_to_cell_glue)
  face_s_nref = lazy_map(constant_field,face_to_nref)
  face_s_nref_minus = face_s_nref

  # Inverse of the Jacobian transpose
  # panel_maps = get_cell_map(cell_grid)
  # cube2panel_map = lazy_map(p-> grid_maps[p], panel_ids)
  # cell_q_x = lazy_map(∘,cube2panel_map,panel_maps)

  _cell_q_x = get_cell_map(cell_grid)
  cell_q_x = lazy_map(∘,cell_geo_map,_cell_q_x)
  cell_q_Jt = lazy_map(∇,cell_q_x)
  cell_q_invJt = lazy_map(Operation(pinvJt),cell_q_Jt)
  face_q_invJt = lazy_map(Reindex(cell_q_invJt),face_to_cell_glue.face_to_cell)

  # Change of domain
  face_s_q = face_glue.tface_to_mface_map
  face_s_invJt = lazy_map(∘,face_q_invJt,face_s_q)
  face_s_n = lazy_map(Broadcasting(Operation(push_normal)),face_s_invJt,face_s_nref)
  cell_vectors =  Fields.MemoArray(face_s_n)

a = lazy_map(evaluate,face_s_nref_plus,Fill(VectorValue(0.,),length(face_s_nref)))
b = lazy_map(evaluate,face_s_nref_minus,Fill(VectorValue(0.,),length(face_s_nref)))

lazy_map(evaluate,face_s_q,Fill([VectorValue(0.,),VectorValue(1.,)],length(face_s_nref)))

face_cmap_minus = lazy_map(Reindex(cell_q_x),face_to_cell_glue.face_to_cell)
face_face_map_minus = lazy_map(∘,face_cmap_minus,face_s_q)
d = lazy_map(evaluate,face_face_map_minus,Fill([VectorValue(0.,),VectorValue(1.,)],length(face_s_nref)))

e = map(.-,c,d)
map(x->vector_length.(x),e)


n_2 = GenericCellField(cell_vectors,trian,ReferenceDomain())

pts = get_cell_points(trian)
(n_1 + n_2)(pts)./1



################################################################################

panel_model = coarse_parametric_model()
panel_model = Adaptivity.refine(panel_model)

btrian = BoundaryTriangulation(panel_model,1)
glue = get_glue(btrian,Val(2))
panel_ids = get_panel_ids(panel_model)

lazy_map(Reindex(panel_ids),glue.tface_to_mface) == panel_ids[glue.tface_to_mface]


############ left and right normals
## left
btrian = BoundaryTriangulation(panel_model,1)
n = get_normal_vector(btrian)
face_panel_ids = get_panel_ids(btrian)

jacobian_cf = panelwise_cellfield(forward_jacobian,btrian,face_panel_ids)
inv_metric_cf = CellField(analytic_inv_metric,btrian)

_n_mapped = jacobian_cf ⋅ (inv_metric_cf ⋅ n)
ff = Operation(sqrt)(  n ⋅ (inv_metric_cf ⋅ n)  )
n_mapped = _n_mapped/ff


panel_cfs = [jacobian_cf, inv_metric_cf, n , n_mapped]
labels = ["jacobian","metric","n", "n_mapped"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)

face_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), face_panel_ids)
writevtk(btrian,dir*"/ambient_model_skeleton_left",cellfields=cellfields,append=false,geo_map=face_geo_map)

## right
btrian = BoundaryTriangulation(panel_model,2)
n = get_normal_vector(btrian)
face_panel_ids = get_panel_ids(btrian)

jacobian_cf = panelwise_cellfield(forward_jacobian,btrian,face_panel_ids)
inv_metric_cf = CellField(analytic_inv_metric,btrian)

_n_mapped = jacobian_cf ⋅ (inv_metric_cf ⋅ n)
ff = Operation(sqrt)(  n ⋅ (inv_metric_cf ⋅ n)  )
n_mapped = _n_mapped/ff


panel_cfs = [jacobian_cf, inv_metric_cf, n , n_mapped]
labels = ["jacobian","metric","n", "n_mapped"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)

face_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), face_panel_ids)
writevtk(btrian,dir*"/ambient_model_skeleton_right",cellfields=cellfields,append=false,geo_map=face_geo_map)


######## skeleton

panel_model = coarse_parametric_model()
# panel_model = Adaptivity.refine(panel_model)

skel_trian = SkeletonTriangulation(panel_model)
n_skel = get_normal_vector(skel_trian)
skel_panel_ids = get_panel_ids(skel_trian)

panel_cfs = [n_skel.plus , n_skel.minus]
labels = ["n_plus", "n_mins"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
skel_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), skel_panel_ids.plus)
writevtk(skel_trian,dir*"/ambient_model_skeleton",cellfields=cellfields,append=false,geo_map=skel_geo_map)


jacobian_cf = panelwise_cellfield(forward_jacobian,skel_trian,skel_panel_ids)
inv_metric_cf = CellField(analytic_inv_metric,skel_trian)

_n_mapped_plus = jacobian_cf.plus ⋅ (inv_metric_cf ⋅ n_skel.plus)
ff = Operation(sqrt)(  n_skel.plus ⋅ (inv_metric_cf ⋅ n_skel.plus)  )
n_mapped_plus = _n_mapped_plus/ff

_n_mapped_minus = jacobian_cf.minus ⋅ (inv_metric_cf ⋅ n_skel.minus)
ff = Operation(sqrt)(  n_skel.minus ⋅ (inv_metric_cf ⋅ n_skel.minus)  )
n_mapped_minus = _n_mapped_minus/ff

pts = get_cell_points(skel_trian)

(n_skel.plus  + n_skel.minus)(pts)

(n_mapped_plus + n_mapped_minus)(pts)




panel_cfs = [n_skel.plus , n_skel.minus, n_skel.minus+n_skel.plus]
labels = ["n_plus", "n_mins","tot"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
skel_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), skel_panel_ids.plus)
writevtk(skel_trian,dir*"/ambient_model_skeleton",cellfields=cellfields,append=false,geo_map=skel_geo_map)




model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1) ,(3,3),isperiodic=(true,true)))
_skel = SkeletonTriangulation(model)
_n = get_normal_vector(_skel)
_pts = get_cell_points(_skel)
(_n.plus + _n.minus)(_pts)


n_skel⋅n_skel
n_skel.plus ⋅ n_skel.plus

cache = nothing
evaluate!(cache,⋅,n_skel.plus,n_skel.plus)
evaluate!(cache,⋅,n_skel,n_skel)



jacobian_cf ⋅n_skel

_n_mapped = jacobian_cf ⋅ (inv_metric_cf ⋅ n_skel)
ff = Operation(sqrt)(  n ⋅ (inv_metric_cf ⋅ n_skel)  )
n_mapped = _n_mapped/ff

############## only 1 panel
skel_trian = SkeletonTriangulation(panel_model)
_skel_panel_ids = get_panel_ids(skel_trian)


face_to_mask = collect(_skel_panel_ids.plus.==1)
skel_trian = SkeletonTriangulation(panel_model,face_to_mask)

n_skel = get_normal_vector(skel_trian)
skel_panel_ids = get_panel_ids(skel_trian)


jacobian_cf = panelwise_cellfield(forward_jacobian,skel_trian,skel_panel_ids)
inv_metric_cf = CellField(analytic_inv_metric,skel_trian)

_n_mapped = jacobian_cf.plus ⋅ (inv_metric_cf ⋅ n_skel.plus)
ff = Operation(sqrt)(  n_skel.plus ⋅ (inv_metric_cf ⋅ n_skel.plus)  )
n_mapped = _n_mapped/ff

_n_mapped = jacobian_cf.minus ⋅ (inv_metric_cf ⋅ n_skel.minus)
ff = Operation(sqrt)(  n_skel.minus ⋅ (inv_metric_cf ⋅ n_skel.minus)  )
n_mapped_minus = _n_mapped/ff

panel_cfs = [n_mapped, n_mapped_minus]
labels = ["n_plus", "n_minus"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
skel_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), skel_panel_ids.plus)
writevtk(skel_trian,dir*"/ambient_model_skeleton",cellfields=cellfields,append=false,geo_map=skel_geo_map)


panel_cfs = [n_skel.plus, n_skel.minus]
labels = ["n_plus", "n_minus"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
skel_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), skel_panel_ids.plus)
writevtk(skel_trian,dir*"/ambient_model_skeleton_flat",cellfields=cellfields,append=false)
