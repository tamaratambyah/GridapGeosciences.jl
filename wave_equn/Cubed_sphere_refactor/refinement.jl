import Gridap.Adaptivity: EdgeBasedRefinement
function Gridap.Adaptivity.refine(model::CubedSphereDiscreteModel,args...;refinement_method="red_green",kwargs...)
  return refine(Gridap.Adaptivity.string_to_refinement(refinement_method, model),model,args...;kwargs...)
end


function Gridap.Adaptivity.refine(method::EdgeBasedRefinement,model::CubedSphereDiscreteModel{Dc,Dp};cells_to_refine=nothing) where {Dc,Dp}

  ctopo = get_grid_topology(model)
  coarse_labels = get_face_labeling(model)
  # Create new model
  rrules, faces_list = Gridap.Adaptivity.setup_edge_based_rrules(method, model.grid_topology,cells_to_refine)
  topo   = Gridap.Adaptivity.refine_edge_based_topology(ctopo,rrules,faces_list)
  reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
  grid   = UnstructuredGrid(
    get_vertex_coordinates(topo),
    get_faces(topo,Dc,0),reffes,
    get_cell_type(topo),
    OrientationStyle(topo)
  )

  cube_to_sphere_map = get_cube_to_sphere_map(get_grid(model))
  if typeof(cube_to_sphere_map) == FEFunction
    refine_cube_to_sphere_map()
  end
  CSgrid = CubedSphereGrid(grid,cube_to_sphere_map)

  glue = Gridap.Adaptivity.blocked_refinement_glue(rrules)
  labels = Gridap.Adaptivity.refine_face_labeling(coarse_labels,glue,model.grid_topology,topo)
  ref_model = CubedSphereDiscreteModel(CSgrid,topo,labels)
  return AdaptedDiscreteModel(ref_model,model,glue)
end
