

panel_model = coarse_parametric_model()
panel_model = Adaptivity.refine(panel_model)



amodel = panel_model

model = get_model(amodel)
# model = panel_model


panel_ids = get_panel_ids(model)
topo = get_grid_topology(amodel)

Dc = num_cell_dims(model)
d = Dc - 1
get_faces(topo,Dc,1)

cmaps = get_cell_map(get_grid(model))

face_to_mask = collect(Bool, .!get_isboundary_face(topo,Dc-1))

# arbitary boundary triangulation to get F2Fglue
# Note, F2Fglue only on reference cell, so okay to use junk nodes
_btrian = BoundaryTriangulation(model,face_to_mask)
F2Fglue = Geometry.get_glue(_btrian,Val(2),Val(2))

face_2_cell = F2Fglue.tface_to_mface

ref_face_2_ref_cell_map = F2Fglue.tface_to_mface_map
cfmaps = map(f-> cmaps[f],face_2_cell)

fmaps = lazy_map(∘, cfmaps,ref_face_2_ref_cell_map)

pts = get_cell_ref_coordinates(trian)

trian_coords = lazy_map(evaluate,fmaps,pts)


node_coordinates = collect1d(get_node_coordinates(model))
cell_to_nodes = Table(get_face_nodes(model,d))
cell_to_type = collect1d(get_face_type(model,d))
reffes = get_reffaces(ReferenceFE{d},model)

bgface_grid = Geometry.UnstructuredGrid(node_coordinates,cell_to_nodes,reffes,cell_to_type,OrientationStyle(cell_grid),
            nothing,fmaps)

face_to_bgface =  findall(face_to_mask)
bgface_to_lcell = fill(1,num_facets(model))

face_grid = view(bgface_grid,face_to_bgface)
cell_grid = get_grid(model)
glue = Geometry.FaceToCellGlue(topo,cell_grid,face_grid,face_to_bgface,bgface_to_lcell)
trian = BodyFittedTriangulation(model,face_grid,face_to_bgface)

bound_trian = BoundaryTriangulation(trian,glue)


atrian = Adaptivity.AdaptedTriangulation(bound_trian,amodel)

Geometry.get_glue(btrian::BoundaryTriangulation) = btrian.glue
Geometry.get_glue(btrian::AdaptedTriangulation) = atrian.trian.glue

glue = get_glue(atrian)
face_2_cell = glue.face_to_cell
face_panel_ids = panel_ids[face_2_cell]
face_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), face_panel_ids)

writevtk(bound_trian,dir*"/ambient_model_skeleton",append=false,geo_map=face_geo_map)





Λ_panel = SkeletonTriangulation(panel_model)
n_Λ = get_normal_vector(Λ_panel)
pts = get_cell_points(Λ_panel)

visdata, = visualization_data(bound_trian,dir*"/ambient_model_skeleton")
get_cell_coordinates(visdata.grid)
