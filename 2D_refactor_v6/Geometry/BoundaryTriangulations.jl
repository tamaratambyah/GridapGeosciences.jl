function Geometry.BoundaryTriangulation(model::ParametricDiscreteModel,lcell::Integer=1)
  topo = get_grid_topology(model)
  Dc = num_cell_dims(model)

  bgface_to_mask = collect(Bool, .!get_isboundary_face(topo,Dc-1))
  bgface_to_lcell = Fill(lcell,num_facets(model))

  Geometry.BoundaryTriangulation(model,bgface_to_mask,bgface_to_lcell)
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

  ## make bgface_grid with proper face maps
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
