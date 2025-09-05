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


panel_model = coarse_parametric_model()
panel_model = Adaptivity.refine(panel_model)
btrian = BoundaryTriangulation(panel_model,2)


face_panel_ids = get_face_panel_ids(btrian)
face_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), face_panel_ids)
writevtk(btrian,dir*"/ambient_model_skeleton",append=false,geo_map=face_geo_map)

function get_face_panel_ids(atrian::AdaptedTriangulation)
  get_face_panel_ids(atrian.trian)
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



### adapted model
panel_model = Adaptivity.refine(panel_model)
panel_ids = get_panel_ids(panel_model)
btrian = BoundaryTriangulation(panel_model,2)

glue = get_glue(btrian,Val(2))
face_2_cell = glue.tface_to_mface
face_panel_ids = panel_ids[face_2_cell]
face_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), face_panel_ids)
writevtk(btrian,dir*"/ambient_model_skeleton",append=false,geo_map=face_geo_map)



### adapted, adapted model
panel_model = Adaptivity.refine(panel_model)
panel_ids = get_panel_ids(panel_model)
btrian = BoundaryTriangulation(panel_model,1)

glue = get_glue(btrian,Val(2))
face_2_cell = glue.tface_to_mface
face_panel_ids = panel_ids[face_2_cell]
face_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), face_panel_ids)

writevtk(btrian,dir*"/ambient_model_skeleton",append=false,geo_map=face_geo_map)



_skel_trian = SkeletonTriangulation(panel_model)

skel_trian = _skel_trian.trian
left_btrian = skel_trian.minus
glue = get_glue(left_btrian,Val(2))
face_2_cell = glue.tface_to_mface
face_panel_ids = panel_ids[face_2_cell]
face_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), face_panel_ids)
writevtk(left_btrian,dir*"/ambient_model_skeleton_minus",append=false,geo_map=face_geo_map)


right_btrian = skel_trian.plus
glue = get_glue(right_btrian,Val(2))
face_2_cell = glue.tface_to_mface
face_panel_ids = panel_ids[face_2_cell]
face_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), face_panel_ids)
writevtk(right_btrian,dir*"/ambient_model_skeleton_plus",append=false,geo_map=face_geo_map)




visdata, = visualization_data(skel_trian,dir*"/ambient_model_skeleton")
get_cell_coordinates(visdata.grid)
