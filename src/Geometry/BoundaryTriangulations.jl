
function Geometry.BoundaryTriangulation(model::CubedSphereParametricDiscreteModel,lcell::Integer=1)
  topo = get_grid_topology(model)
  Dc = num_cell_dims(model)

  bgface_to_mask = collect(Bool, .!get_isboundary_face(topo,Dc-1))
  bgface_to_lcell = Fill(lcell,num_facets(model))

  Geometry.BoundaryTriangulation(model,bgface_to_mask,bgface_to_lcell)
end


function Geometry.BoundaryTriangulation(
  model::CubedSphereParametricDiscreteModel,
  bgface_to_mask::AbstractVector{Bool},
  bgface_to_lcell::AbstractVector{<:Integer}
  )
  face_to_bgface =  findall(bgface_to_mask)
  Geometry.BoundaryTriangulation(model,face_to_bgface,bgface_to_lcell)
end



function Geometry.BoundaryTriangulation(
  model::CubedSphereParametricDiscreteModel,
  face_to_bgface::AbstractVector{<:Integer},
  bgface_to_lcell::AbstractVector{<:Integer}
  )

  # println("my boundary trian")
  # println(bgface_to_lcell)

  topo = get_grid_topology(model)
  Dc = num_cell_dims(model)
  d = Dc - 1

  # Note, restrict the cell_to_nodes and cell_to_type based on the mask
  # This is because the Distributed constructor dispatches here and in Distributed,
  # we do not always have a left cell / right cell on every processor
  node_coordinates = collect1d(get_node_coordinates(model))
  cell_to_nodes = Table(lazy_map(Reindex(get_face_nodes(model,d)),face_to_bgface))
  cell_to_type = collect1d(lazy_map(Reindex(get_face_type(model,d)),face_to_bgface))
  reffes = get_reffaces(ReferenceFE{d},model)

  # make the facegrid with wrong fmaps to get F2Fglue
  # Note, F2Fglue only on reference cell, so okay to use junk nodes
  _face_grid = Geometry.UnstructuredGrid(node_coordinates,cell_to_nodes,reffes,cell_to_type,NonOriented())
  _cell_grid = get_grid(model)
  _btrian = BodyFittedTriangulation(model,_face_grid,face_to_bgface)
  _glue = Geometry.FaceToCellGlue(topo,_cell_grid,_face_grid,face_to_bgface,bgface_to_lcell)
  _trian = BoundaryTriangulation(_btrian,_glue)
  F2Fglue = Geometry.get_glue(_trian,Val(Dc),Val(Dc))
  # println(Dc)
  # tface_to_mface_map = Gridap.Geometry.compute_face_to_cell_reference_map(_cell_grid,_face_grid,_glue)

  ## get the face2cell_ids
  face_2_cell_ids = _glue.face_to_cell; #println(length(face_2_cell_ids))

  # compose face map and cell map
  ref_face_2_ref_cell_map = F2Fglue.tface_to_mface_map;   #println((ref_face_2_ref_cell_map))
  cmaps = get_cell_map(get_grid(model))

  # cfmaps = map(f-> cmaps[f],face_2_cell_ids)
  cfmaps = lazy_map(Reindex(cmaps),face_2_cell_ids)
  fmaps = lazy_map(∘, cfmaps,ref_face_2_ref_cell_map)

  # make the facegrid with correct fmaps
  face_grid = Geometry.UnstructuredGrid(node_coordinates,cell_to_nodes,reffes,cell_to_type,NonOriented(),nothing,fmaps)
  cell_grid = get_grid(model)
  glue = Geometry.FaceToCellGlue(topo,cell_grid,face_grid,face_to_bgface,bgface_to_lcell)
  trian = BodyFittedTriangulation(model,face_grid,face_to_bgface)

  Geometry.BoundaryTriangulation(trian,glue)

end


################################################################################
"""
get the facet normal vector on the sphere
compute as the pushforward of the parametric facet normal
Two different methods:
  1. using Gridap's machinary
  2. using Santi's formula
"""

function pushforward_normal(args...)
  @abstractmethod
end

"""
map normal vector using cell_geo_map
This method is an adaption of Gridap's machinary
"""
function pushforward_normal(trian::BoundaryTriangulation,cell_geo_map::AbstractArray)
  cell_vectors = get_facet_normal(trian,cell_geo_map)
  get_normal_vector(trian,cell_vectors)
end

Geometry.get_facet_normal(trian::BoundaryTriangulation,cell_geo_map::AbstractArray) = Geometry.get_facet_normal(trian,trian.glue,cell_geo_map)

function Geometry.get_facet_normal(trian::BoundaryTriangulation,face_to_cell_glue,cell_geo_map::AbstractArray)
  bgmodel = get_background_model(trian)
  D = num_cell_dims(bgmodel)
  cell_grid = get_grid(bgmodel)
  face_glue = get_glue(trian,Val(D))

  @check length(cell_geo_map) == num_cells(cell_grid) "\n cell_geo_map must be a cell-wise array"

  ## Map to ambient space
  return get_mapped_facet_normal(cell_grid,face_glue,face_to_cell_glue,cell_geo_map)

end

function get_mapped_facet_normal(
  cell_grid::Grid,
  face_glue::FaceToFaceGlue,
  face_to_cell_glue::FaceToCellGlue,
  cell_geo_map::AbstractArray
)
  # println("mapped normal")

  @check length(cell_geo_map) == num_cells(cell_grid) "\n cell_geo_map must be a cell-wise array
  "

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
  Gridap.Fields.MemoArray(face_s_n)
end



"""
push normal vector from chart to the surface
This method is based on Santi's formula
"""
pushforward_normal(trian::BoundaryTriangulation) = _pushforward_normal(trian)
function _pushforward_normal(trian)
  panel_model = get_background_model(trian)
  panel_model_metadata = get_forward_map_generator(panel_model)
  n_2_2D = get_normal_vector(trian)

  face_panel_ids = get_panel_ids(trian)
  glue = get_glue(trian,Val(2))

  face_geo_map = lazy_map(p -> panel_model_metadata(p), face_panel_ids)
  Jt = lazy_map(∇,face_geo_map)
  J = lazy_map(Operation(transpose),Jt)

  cell_q_x = get_cell_map(get_grid(panel_model))
  _cell_q_x = lazy_map(Reindex(cell_q_x),glue.tface_to_mface)
  cell_pts = lazy_map(∘,_cell_q_x,glue.tface_to_mface_map)
  J_face = lazy_map(∘,J,cell_pts)

  J_cf = GenericCellField(Gridap.Fields.MemoArray(J_face),trian,ReferenceDomain())

  # inv_cf = CellField(analytic_inv_metric,trian)
  inv_cf = ParametricCellField(inv_metric,trian)

  _n_mapped = J_cf ⋅ (inv_cf  ⋅ n_2_2D )
  ff = Operation(sqrt)(  n_2_2D   ⋅ (inv_cf⋅ n_2_2D )  )
  n_mapped = _n_mapped/ff

  return n_mapped
end


"""
pullback of the area form
returns |Jg^{-1} ̂n|
"""
##
pullback_area_form(trian::BoundaryTriangulation) = _pullback_area_form(trian)
function _pullback_area_form(trian)
  inv_metric_cf = ParametricCellField(inv_metric,trian)
  jac_cf = ParametricCellField(forward_jacobian,trian)
  n_Λ = get_normal_vector(trian)

  Operation(norm)(jac_cf⋅(inv_metric_cf ⋅n_Λ) )
end



"""
BoundaryTriangulation

For the CubedSphereAmbientDiscreteModel, the boundary triangulation is the pushforward of the
boundary triangulation of the underlying parametric model. To achieve this,
  1. create boundary triangulation of the parametric model
  2. create a ambient_face_grid that is the composition of the panel_cmap and the
  forward map
"""
function Geometry.BoundaryTriangulation(a_model::CubedSphereAmbientDiscreteModel;tags=nothing)
  labeling = get_face_labeling(a_model)
  model = get_parametric_model(a_model)
  btrian = BoundaryTriangulation(model,labeling,tags=tags)
  pushforward_trian(a_model,btrian)
end

function Geometry.BoundaryTriangulation(
  a_model::CubedSphereAmbientDiscreteModel, bgface_to_mask::AbstractVector{Bool}, lcell::Integer=1)
  model = get_parametric_model(a_model)
  btrian = BoundaryTriangulation(model,bgface_to_mask, Fill(lcell,num_facets(model)) )
  pushforward_trian(a_model,btrian)
end

function Geometry.SkeletonTriangulation(a_model::CubedSphereAmbientDiscreteModel;tags=nothing)
  labeling = get_face_labeling(a_model)
  model = get_parametric_model(a_model)
  strian = SkeletonTriangulation(model,labeling;tags=tags)
  plus = pushforward_trian(a_model,strian.plus)
  minus = pushforward_trian(a_model,strian.minus)
  SkeletonTriangulation(plus,minus)
end

function pushforward_trian(model::CubedSphereAmbientDiscreteModel,panel_trian)

  panel_model = get_background_model(panel_trian)

  # Extract the panel_face_grid and cmaps that go to alpha,beta.
  # compose alpha,beta -> manifold to get ambient_fmap
  panel_face_grid = get_grid(panel_trian)
  panel_cmap = get_cell_map(panel_face_grid)

  fwd_map_generator = get_forward_map_generator(panel_model)
  panel_ids = get_panel_ids(panel_trian)
  fwd_map =  geo_map_func(fwd_map_generator,panel_ids)

  ambient_cmap = lazy_map(∘,fwd_map,panel_cmap)

  # Create the ambient_face_grid with the correct cmap
  D = num_cell_dims(panel_model)
  d = D - 1

  node_coordinates = Gridap.Arrays.collect1d(Gridap.Geometry.get_node_coordinates(model))
  cell_to_nodes = get_cell_node_ids(panel_trian)
  cell_to_type = Gridap.Geometry.get_cell_type(panel_face_grid)
  reffes = Gridap.Geometry.get_reffaces(ReferenceFE{d},model)

  ambient_face_grid = Geometry.UnstructuredGrid(node_coordinates,cell_to_nodes,reffes,cell_to_type,NonOriented(),nothing,ambient_cmap)

  glue = panel_trian.glue
  face_to_bgface = panel_trian.glue.face_to_bgface
  trian = BodyFittedTriangulation(model,ambient_face_grid,face_to_bgface)

  Geometry.BoundaryTriangulation(trian,glue)
end
