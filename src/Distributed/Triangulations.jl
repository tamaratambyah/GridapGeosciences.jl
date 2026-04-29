# function GridapDistributed.BoundaryTriangulation(
#   portion,model::CubedSphere2DParametricDistributedDiscreteModel,labels::DistributedFaceLabeling;tags=nothing)
#   println("distributed booundary trian")
#   Dc = num_cell_dims(model)

#   topo = get_grid_topology(model)
#   face_to_mask = get_isboundary_face(topo,Dc-1) # This is globally consistent

#   ## for ParametricDiscreteModel, we want all cells all the internal cells
#   _face_to_mask = map(face_to_mask) do m
#     return .!m
#   end

#   Gridap.Geometry.BoundaryTriangulation(portion,model,_face_to_mask)
# end


function pullback_area_form(trian::DistributedTriangulation)
  fields = map(trian.trians) do t
    return pullback_area_form(t)
  end
  DistributedCellField(fields,trian)
end


function pushforward_normal(trian::GridapDistributed.DistributedTriangulation,cell_geo_map::AbstractArray)
  fields = map(trian.trians,cell_geo_map) do t,m
    pushforward_normal(t,m)
  end
  GridapDistributed.DistributedCellField(fields,trian)
end

function pushforward_normal(trian::GridapDistributed.DistributedTriangulation)
  fields = map(trian.trians) do t
    return pushforward_normal(t)
  end
  return GridapDistributed.DistributedCellField(fields,trian)
end
