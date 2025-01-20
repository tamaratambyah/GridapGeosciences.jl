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
  cube_grid   = UnstructuredGrid(
    get_vertex_coordinates(topo),
    get_faces(topo,Dc,0),reffes,
    get_cell_type(topo),
    OrientationStyle(topo)
  )

  glue = Gridap.Adaptivity.blocked_refinement_glue(rrules)
  labels = Gridap.Adaptivity.refine_face_labeling(coarse_labels,glue,model.grid_topology,topo)
  _ref_model = UnstructuredDiscreteModel(cube_grid,topo,labels)
  cube_modelh = AdaptedDiscreteModel(_ref_model,model,glue)


  cube_to_sphere_map = get_cube_to_sphere_map(get_grid(model))
  println( typeof(cube_to_sphere_map) <: FEFunction )
  if typeof(cube_to_sphere_map) <: FEFunction
    maph,order,Vh = transfer_FE_map(cube_modelh,cube_to_sphere_map)
    return CSgrid = CubedSphereGrid(cube_modelh,maph,order)
  elseif typeof(cube_to_sphere_map) <: Function
    return CSgrid = CubedSphereGrid(cube_grid,cube_to_sphere_map)
  end

  ref_model = CubedSphereDiscreteModel(CSgrid,topo,labels)
  return AdaptedDiscreteModel(ref_model,model,glue)
end
