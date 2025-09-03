
function BoundaryTriangulation(model::ParametricDiscreteModel,lcell::Int=1)

  topo = get_grid_topology(model)
  Dc = num_cell_dims(model)
  d = Dc - 1

  # arbitary boundary triangulation to get F2Fglue
  # Note, F2Fglue only on reference cell, so okay to use junk nodes
  face_to_mask = collect(Bool, .!get_isboundary_face(topo,Dc-1))
  _btrian = BoundaryTriangulation(model,face_to_mask)
  F2Fglue = Geometry.get_glue(_btrian,Val(2),Val(2))

  ## get the face2cell_ids
  face_2_cell_ids = F2Fglue.tface_to_mface

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
  face_to_bgface =  findall(face_to_mask)
  bgface_to_lcell = fill(lcell,num_facets(model))

  face_grid = view(bgface_grid,face_to_bgface)
  cell_grid = get_grid(model)
  glue = Geometry.FaceToCellGlue(topo,cell_grid,face_grid,face_to_bgface,bgface_to_lcell)
  trian = BodyFittedTriangulation(model,face_grid,face_to_bgface)

  BoundaryTriangulation(trian,glue)

end

# panel_model = Adaptivity.refine(panel_model)
# amodel = panel_model

# model = get_model(amodel)
# model = panel_model

# get_faces(topo,Dc,1)

panel_model = coarse_parametric_model()
panel_ids = get_panel_ids(panel_model)
btrian = BoundaryTriangulation(panel_model,2)

glue = get_glue(btrian,Val(2))
face_2_cell = glue.tface_to_mface
face_panel_ids = panel_ids[face_2_cell]
face_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), face_panel_ids)
writevtk(btrian,dir*"/ambient_model_skeleton",append=false,geo_map=face_geo_map)

### adapted model
panel_model = Adaptivity.refine(panel_model)
panel_ids = get_panel_ids(panel_model)
btrian = BoundaryTriangulation(panel_model,1)

glue = get_glue(btrian,Val(2))
face_2_cell = glue.tface_to_mface
face_panel_ids = panel_ids[face_2_cell]
face_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), face_panel_ids)
writevtk(btrian,dir*"/ambient_model_skeleton",append=false,geo_map=face_geo_map)



### adapted, adapted model
panel_model = Adaptivity.refine(panel_model)
panel_ids = get_panel_ids(panel_model)
btrian = BoundaryTriangulation(panel_model,2)

glue = get_glue(btrian,Val(2))
face_2_cell = glue.tface_to_mface
face_panel_ids = panel_ids[face_2_cell]
face_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), face_panel_ids)

writevtk(btrian,dir*"/ambient_model_skeleton",append=false,geo_map=face_geo_map)



function SkeletonTriangulation(model::ParametricDiscreteModel)
  println("my skeleton")
  left_cell_around = 1
  plus = BoundaryTriangulation(model,left_cell_around)
  right_cell_around = 2
  minus = BoundaryTriangulation(model,right_cell_around)
  SkeletonTriangulation(plus,minus)
end


skel_trian = SkeletonTriangulation(panel_model)

writevtk(skel_trian,dir*"/ambient_model_skeleton",append=false,geo_map=face_geo_map)


visdata, = visualization_data(skel_trian,dir*"/ambient_model_skeleton")
get_cell_coordinates(visdata.grid)
