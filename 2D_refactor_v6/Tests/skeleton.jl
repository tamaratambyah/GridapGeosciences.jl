panel_model = coarse_parametric_model()
# panel_model = Gridap.Adaptivity.refine(panel_model)

panel_ids = get_panel_ids(panel_model)






topo = get_grid_topology(panel_model)

D = 2
get_faces(topo,D,0)
get_faces(topo,D,1)
get_faces(topo,D,2)





Λ_panel = SkeletonTriangulation(panel_model)
n_Λ = get_normal_vector(Λ_panel)
pts = get_cell_points(Λ_panel)

visdata, = visualization_data(Λ_panel,dir*"/ambient_model_skeleton")
get_cell_coordinates(visdata.grid)


skeleton_panel_ids = [1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5]
skeleton_cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), skeleton_panel_ids)


skel_ref_points = get_cell_ref_coordinates(Λ_panel)
skel_edges =  lazy_map(evaluate,skeleton_cell_geo_map,)


writevtk(Λ_panel,dir*"/ambient_model_skeleton",append=false,geo_map=skeleton_cell_geo_map)

jacobian_cf = panelwise_cellfield(forward_jacobian,Λ_panel,panel_ids)

physical_normals = jacobian_cf.minus ⋅ n_Λ.minus

writevtk(Λ_panel,dir*"/ambient_model_skeleton",cellfields=["nplus"=>physical_normals],append=false)


topo = get_grid_topology(panel_model)
D = num_cell_dims(panel_model)
face_to_mask = collect(Bool, .!get_isboundary_face(topo,D-1))




########## test 2D model
nodes  =   [
    Point(0.0, 0.0)  # node 1
    Point(1.0, 0.0)   # node 2
    Point(0.0, 1.0)   # node 3
    Point(1.0, 1.0)    # node 4
    Point(2.0, 0.0)  # node 5
    Point(2.0, 1.0)   # node 6

  ]

  ## CCAM panel ordering
  data = [ 1,2,3,4, 2,5,4,6 ]

  ptr = generate_ptr(2)
  cell_node_ids = Table(data,ptr)

  polytopes = fill(QUAD,2)
  cell_type = fill(1,2)
  reffes = LagrangianRefFE(Float64,QUAD,1)
  cell_reffes=[reffes]

  topo = UnstructuredGridTopology(nodes,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
  labels = FaceLabeling(topo)

  cube_grid = Gridap.Geometry.UnstructuredGrid(nodes,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

model = UnstructuredDiscreteModel(cube_grid,topo,labels)
Λ = SkeletonTriangulation(model)
writevtk(Λ,dir*"/cube",append=false)


topo = get_grid_topology(model)
D = num_cell_dims(model)
face_to_mask = collect(Bool, .!get_isboundary_face(topo,D-1))
