function Geometry.BoundaryTriangulation(model::ParametricDiscreteModel,lcell=1)
  topo = get_grid_topology(model)
  Dc = num_cell_dims(model)

  face_to_mask = collect(Bool, .!get_isboundary_face(topo,Dc-1))
  face_to_bgface =  findall(face_to_mask)
  bgface_to_lcell = Fill(lcell,num_facets(model))

  Geometry.BoundaryTriangulation(model,face_to_bgface,bgface_to_lcell)
end


function Geometry.BoundaryTriangulation(
  model::ParametricDiscreteModel,
  bgface_to_mask::AbstractVector{Bool},
  bgface_to_lcell::AbstractVector{<:Integer}
  )
  face_to_bgface =  findall(bgface_to_mask)
  Geometry.BoundaryTriangulation(model,face_to_bgface,bgface_to_lcell)
end

function Geometry.BoundaryTriangulation(
  model::ParametricDiscreteModel,
  face_to_bgface::AbstractVector{<:Integer},
  bgface_to_lcell::AbstractVector{<:Integer}
  )
  println("my boundary trian")

  topo = get_grid_topology(model)
  Dc = num_cell_dims(model)
  d = Dc - 1

  # arbitary boundary triangulation to get F2Fglue
  # Note, F2Fglue only on reference cell, so okay to use junk nodes
  _bgface_grid = Grid(ReferenceFE{Dc-1},model)
  _face_grid = view(_bgface_grid,face_to_bgface)
  _cell_grid = get_grid(model)
  _btrian = BodyFittedTriangulation(model,_face_grid,face_to_bgface)
  _glue = Geometry.FaceToCellGlue(topo,_cell_grid,_face_grid,face_to_bgface,bgface_to_lcell)
  _trian = BoundaryTriangulation(_btrian,_glue)
  F2Fglue = Geometry.get_glue(_trian,Val(2),Val(2))

  # tface_to_mface_map = Gridap.Geometry.compute_face_to_cell_reference_map(_cell_grid,_face_grid,_glue)

  ## get the face2cell_ids
  face_2_cell_ids = _glue.face_to_cell

  # compose face map and cell map
  ref_face_2_ref_cell_map = F2Fglue.tface_to_mface_map
  cmaps = get_cell_map(get_grid(model))
  cfmaps = map(f-> cmaps[f],face_2_cell_ids)
  fmaps = lazy_map(∘, cfmaps,ref_face_2_ref_cell_map)

  # pts = get_cell_ref_coordinates(trian)
  # trian_coords = lazy_map(evaluate,fmaps,pts)

  ## make bgface_grid with proper cells maps
  node_coordinates = collect1d(get_node_coordinates(model))
  cell_to_nodes = Table(get_face_nodes(model,d))
  cell_to_type = collect1d(get_face_type(model,d))
  reffes = get_reffaces(ReferenceFE{d},model)
  bgface_grid = Geometry.UnstructuredGrid(node_coordinates,cell_to_nodes,reffes,cell_to_type,NonOriented(),
              nothing,fmaps)

  # make the boundary triangulation
  face_grid = view(bgface_grid,face_to_bgface)
  cell_grid = get_grid(model)
  glue = Geometry.FaceToCellGlue(topo,cell_grid,face_grid,face_to_bgface,bgface_to_lcell)
  trian = BodyFittedTriangulation(model,face_grid,face_to_bgface)

  Geometry.BoundaryTriangulation(trian,glue)

end

function get_panel_ids(atrian::AdaptedTriangulation)
  get_panel_ids(atrian.trian)
end

function get_panel_ids(btrian::BoundaryTriangulation)
  get_face_panel_ids(btrian)
end

function get_panel_ids(strian::SkeletonTriangulation)
  plus = get_face_panel_ids(strian.plus)
  minus = get_face_panel_ids(strian.plus)
  SkeletonPair(plus,minus)
end

function get_face_panel_ids(btrian::BoundaryTriangulation)
  panel_model = get_background_model(btrian)
  panel_ids = get_panel_ids(panel_model)
  Dc = num_cell_dims(panel_model)

  glue = get_glue(btrian,Val(Dc))
  face_2_cell = glue.tface_to_mface
  face_panel_ids = panel_ids[face_2_cell]
  return face_panel_ids
end





panel_model = coarse_parametric_model()
panel_model = Adaptivity.refine(panel_model)

btrian = BoundaryTriangulation(panel_model,1)

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
panel_model = Adaptivity.refine(panel_model)

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

inv_metric_cf ⋅n_skel

skel_trian = SkeletonTriangulation(panel_model)
n_skel = get_normal_vector(skel_trian)
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
