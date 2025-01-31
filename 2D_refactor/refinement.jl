import Gridap.Adaptivity: EdgeBasedRefinement
function Gridap.Adaptivity.refine(model::ManifoldDiscreteModel,args...;refinement_method="red_green",kwargs...)
  return refine(Gridap.Adaptivity.string_to_refinement(refinement_method, model),model,args...;kwargs...)
end


function Gridap.Adaptivity.refine(method::EdgeBasedRefinement,model::ManifoldDiscreteModel{Dc,Dp_grid,Dp_topo};cells_to_refine=nothing) where {Dc,Dp_grid,Dp_topo}

  ctopo = get_grid_topology(model)
  coarse_labels = get_face_labeling(model)
  # Create new model
  rrules, faces_list = Gridap.Adaptivity.setup_edge_based_rrules(method, ctopo, cells_to_refine)
  topo   = Gridap.Adaptivity.refine_edge_based_topology(ctopo,rrules,faces_list)
  reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
  cube_grid   = UnstructuredGrid(
    get_vertex_coordinates(topo),
    get_faces(topo,Dc,0),reffes,
    get_cell_type(topo),
    OrientationStyle(topo)
  )

  glue = Gridap.Adaptivity.blocked_refinement_glue(rrules)
  labels = Gridap.Adaptivity.refine_face_labeling(coarse_labels,glue,ctopo,topo)
  cube_modelh = UnstructuredManifoldDiscreteModel(cube_grid,topo,labels)

  # bump refined cube to sphere
  cube_to_sphere_map = get_cube_to_sphere_map(get_grid(model))
  transfer_info = get_transfer_info(get_grid(model))

  if typeof(cube_to_sphere_map) <: FEFunction
    analytical_cube_to_sphere_map,transfer,order = transfer_info
    # maph,Vh = transfer_FE_map(cube_modelh,cube_to_sphere_map,order)
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

  ref_model = UnstructuredManifoldDiscreteModel(CSgrid,topo,labels)
  return AdaptedManifoldDiscreteModel(ref_model,model,glue)
end


function transfer_FE_map(cube_modelh::ManifoldDiscreteModel,mapH::FEFunction,order::Integer)

  println(order)

  T_vec = eltype(get_node_coordinates(cube_modelh))
  Vh = FESpace(cube_modelh,
              ReferenceFE(lagrangian,T_vec,order),
              conformity=:H1)

  maph = interpolate(mapH,Vh)

  return maph, Vh

end
