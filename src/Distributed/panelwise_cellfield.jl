# function panelwise_cellfield(f::Function,
#   trian::GridapDistributed.DistributedTriangulation,
#   panel_ids::AbstractArray)

#   dmodel = get_background_model(trian)
#   gids = get_cell_gids(dmodel)

#   fields = map(trian.trians,partition(gids),panel_ids) do t,cid,pid
#     owned_cells = own_to_local(cid)
#     panel_ids = pid[owned_cells]
#     panelwise_cellfield(f,t,panel_ids)
#   end
#   GridapDistributed.DistributedCellField(fields,trian)
# end

#### This function is Alberto's solution to the issue with ghost+local on
#### octree periodic meshes. See issue %%%% in GridapP4est
function panelwise_cellfield(f::Function,
  trian::GridapDistributed.DistributedTriangulation,
  panel_ids::AbstractArray)

  # println("new panelwise cellfield")

  panel_model = trian.model

  fields = map(panel_ids,local_views(panel_model)) do panel_ids, lmodel
    panelwise_cellfield(f,Triangulation(lmodel),panel_ids)
  end

  trians = map(local_views(panel_model)) do lmodel
    Triangulation(lmodel)
  end

  _trian = GridapDistributed.DistributedTriangulation(trians,panel_model)
  GridapDistributed.DistributedCellField(fields,_trian)
end


function panelwise_cellfield(f::Function,
  trian::GridapDistributed.DistributedTriangulation)
  # println("old panelwise cellfield -- skeleton")
  fields = map(trian.trians) do t
    panelwise_cellfield(f,t)
  end
  GridapDistributed.DistributedCellField(fields,trian)
end

### This function is Alberto's solution to the issue with ghost+local on
### octree periodic meshes. See issue %%%% in GridapP4est
# function panelwise_cellfield(f::Function,
#   trian::GridapDistributed.DistributedTriangulation)

#   println("new panelwise cellfield for skeleton mesh")

#   panel_model = trian.model

#   fields = map(local_views(panel_model)) do lmodel
#     panelwise_cellfield(f,SkeletonTriangulation(lmodel))
#     # panelwise_cellfield(f,trian)
#   end

#   # trians = map(local_views(panel_model)) do lmodel
#   #   SkeletonTriangulation(lmodel)
#   # end

#   # _trian = GridapDistributed.DistributedTriangulation(trians,panel_model)
#   GridapDistributed.DistributedCellField(fields,trian)
# end


function geo_map_func(trian::DistributedTriangulation)
  # println("distributed geo map")
  model = get_background_model(trian)
  model_metadata = get_forward_map_generator(model)
  owned_panel_ids = get_owned_panel_ids(model)
  return geo_map_func(model_metadata, owned_panel_ids)
end


function geo_map_func(model_metadata,owned_panel_ids::AbstractArray)
  # println("distributed geo map")

  @assert typeof(owned_panel_ids) <: DebugArray || typeof(owned_panel_ids) <: MPIArray "\n Not distributed panel ids"

  cell_geo_map = map(owned_panel_ids) do pid
    return lazy_map(p -> model_metadata(p), pid)
  end
  return cell_geo_map
end

### latlong geo map func
function latlon_geo_map_func(trian::GridapDistributed.DistributedTriangulation)
  # println("distributed trian latlon geo map")

  model = get_background_model(trian)
  model_metadata = get_forward_map_generator(model)
  owned_panel_ids = get_owned_panel_ids(model)
  return latlon_geo_map_func(model_metadata, owned_panel_ids)
end

function latlon_geo_map_func(model_metadata, owned_panel_ids::AbstractArray)
  # println("distributed latlon geo map")

  @assert typeof(owned_panel_ids) <: DebugArray || typeof(owned_panel_ids) <: MPIArray "\n Not distributed panel ids"

  latlon_cell_geo_map = map(owned_panel_ids) do pid
    cell_geo_map = lazy_map(p -> model_metadata(p), pid)
    fi = lazy_map(p->Cartesian2SphereicalMap(),pid)
    return lazy_map(∘, fi, cell_geo_map)
  end
  return latlon_cell_geo_map
end
