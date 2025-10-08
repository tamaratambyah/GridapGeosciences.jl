function GridapDistributed.BoundaryTriangulation(
  portion,model::DistributedParametricDiscreteModel,labels::GridapDistributed.DistributedFaceLabeling;tags=nothing)
  println("distribue booundary trian")
  Dc = num_cell_dims(model)

  topo = get_grid_topology(model)
  face_to_mask = get_isboundary_face(topo,Dc-1) # This is globally consistent

  ## for ParametricDiscreteModel, we want all cells all the internal cells
  _face_to_mask = map(face_to_mask) do m
    return .!m
  end

  Gridap.Geometry.BoundaryTriangulation(portion,model,_face_to_mask)
end
