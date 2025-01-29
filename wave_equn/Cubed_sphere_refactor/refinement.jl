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
  cube_modelh = CubedSphereAdaptedDiscreteModel(_ref_model,model,glue)

  # bump refined cube to sphere
  cube_to_sphere_map = get_cube_to_sphere_map(get_grid(model))
  transfer_info = get_transfer_info(get_grid(model))

  if typeof(cube_to_sphere_map) <: FEFunction
    analytical_cube_to_sphere_map,transfer,order = transfer_info
    maph,Vh = transfer_FE_map(cube_modelh,cube_to_sphere_map,order)
    if transfer # transfer the map
      println("transferring map")
      CSgrid = CubedSphereGrid(cube_modelh,maph,order,transfer_info)
    else # re-inteprolate the map
      println("re-interpolating map")
      CSgrid = CubedSphereGrid(cube_grid,analytical_cube_to_sphere_map,order;transfer=transfer)
    end

  elseif typeof(cube_to_sphere_map) <: Function
    println("analytical map")
    CSgrid = CubedSphereGrid(cube_grid,cube_to_sphere_map)
  end

  ref_model = CubedSphereDiscreteModel(CSgrid,topo,labels)
  return CubedSphereAdaptedDiscreteModel(ref_model,model,glue)
end


function transfer_FE_map(cube_modelh::CubedSphereAdaptedDiscreteModel,mapH::FEFunction,order::Integer)

  println(order)

  T_vec = eltype(get_node_coordinates(cube_modelh))
  Vh = FESpace(cube_modelh,
              ReferenceFE(lagrangian,T_vec,order),
              conformity=:H1)

  maph = interpolate(mapH,Vh)

  return maph, Vh

end
