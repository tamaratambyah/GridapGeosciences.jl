import Gridap.Adaptivity: EdgeBasedRefinement

function Gridap.Adaptivity.refine(method::EdgeBasedRefinement,model::ManifoldDiscreteModel;cells_to_refine=nothing)
  println("RG refinement")
  Dc = num_cell_dims(model)

  # cells_to_refine can be
  #    a) nothing -> All cells get refined
  #    b) AbstractArray{<:Bool} of size num_cells(model)
  #            -> Only cells such that cells_to_refine[iC] == true get refined
  #    c) AbstractArray{<:Integer}
  #            -> Cells for which gid ∈ cells_to_refine get refined

  manifold_grid = get_grid(model)
  cube_grid_3D = get_3D_cube_grid(manifold_grid)

  ctopo = UnstructuredGridTopology(cube_grid_3D) #get_grid_topology(model)
  coarse_labels = get_face_labeling(model)

  # Create new model
  rrules, faces_list = Adaptivity.setup_edge_based_rrules(method, model.grid_topology,cells_to_refine)
  topo   = Adaptivity.refine_edge_based_topology(ctopo,rrules,faces_list)
  reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
  ref_cube_grid_3D   = UnstructuredGrid(
    get_vertex_coordinates(topo),
    get_faces(topo,Dc,0),reffes,
    get_cell_type(topo),
    OrientationStyle(topo)
  )
  glue = Adaptivity.blocked_refinement_glue(rrules)
  labels = Adaptivity.refine_face_labeling(coarse_labels,glue,model.grid_topology,topo)

  # make the manifold grid
  panel_ids = get_panel_ids(model)
  new_to_old_cells = glue.n2o_faces_map[Dc+1]
  ref_panel_ids = panel_ids[new_to_old_cells]
  manifold_grid = ManifoldGrid(get_manifold_name(model),ref_cube_grid_3D,topo,ref_panel_ids)

  ref_model = ManifoldDiscreteModel(manifold_grid,topo,labels)
  adapted_model = AdaptedDiscreteModel(ref_model,model,glue)

  return adapted_model
end

function Adaptivity.refine(model::ManifoldDiscreteModel,args...;refinement_method="red_green",kwargs...)
  return refine(Adaptivity.string_to_refinement(refinement_method, model),model,args...;kwargs...)
end
