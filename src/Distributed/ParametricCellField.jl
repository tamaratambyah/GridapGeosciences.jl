"""
ParametricCellField

Is the distributed implementation of ParametricCellField.
In such function, we call ParametricCellField on the local model and then
recompute the triangulation to ensure proper handling of ghost cells in octree periodic meshes.

"""

function ParametricCellField(f::Function,
  trian::GridapDistributed.DistributedTriangulation)

  panel_model = trian.model

  fields = map(local_views(panel_model)) do lmodel
    ParametricCellField(f,Triangulation(lmodel))
  end

  trians = map(local_views(panel_model)) do lmodel
    Triangulation(lmodel)
  end

  _trian = GridapDistributed.DistributedTriangulation(trians,panel_model)
  GridapDistributed.DistributedCellField(fields,_trian)
end



function geo_map_func(trian::DistributedTriangulation)
  model = get_background_model(trian)
  fwd_map_generator = get_forward_map_generator(model)
  owned_panel_ids = get_owned_panel_ids(model)
  return geo_map_func(fwd_map_generator, owned_panel_ids)
end


function geo_map_func(fwd_map_generator::AbstractArray,owned_panel_ids::AbstractArray)
  @assert typeof(owned_panel_ids) <: DebugArray || typeof(owned_panel_ids) <: MPIArray "\n Not distributed panel ids"
  cell_geo_map = map(owned_panel_ids,fwd_map_generator) do pid, generator
    return lazy_map(p -> generator(p), pid)
  end
  return cell_geo_map
end

### latlong geo map func
function latlon_geo_map_func(trian::GridapDistributed.DistributedTriangulation)
  model = get_background_model(trian)
  fwd_map_generator = get_forward_map_generator(model)
  owned_panel_ids = get_owned_panel_ids(model)
  return latlon_geo_map_func(fwd_map_generator, owned_panel_ids)
end

function latlon_geo_map_func(fwd_map_generator::AbstractArray, owned_panel_ids::AbstractArray)
  @assert typeof(owned_panel_ids) <: DebugArray || typeof(owned_panel_ids) <: MPIArray "\n Not distributed panel ids"
  latlon_cell_geo_map = map(owned_panel_ids,fwd_map_generator) do pid, generator
    cell_geo_map = lazy_map(p -> generator(p), pid)
    fi = lazy_map(p->Cartesian2SphericalMap(),pid)
    return lazy_map(∘, fi, cell_geo_map)
  end
  return latlon_cell_geo_map
end
