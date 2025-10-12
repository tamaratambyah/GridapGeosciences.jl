function panelwise_cellfield(f::Function,
  trian::GridapDistributed.DistributedTriangulation,
  panel_ids::AbstractArray)

  dmodel = get_background_model(trian)
  gids = get_cell_gids(dmodel)

  fields = map(trian.trians,partition(gids),panel_ids) do t,cid,pid
    owned_cells = own_to_local(cid)
    panel_ids = pid[owned_cells]
    panelwise_cellfield(f,t,panel_ids)
  end
  GridapDistributed.DistributedCellField(fields,trian)
end

function panelwise_cellfield(f::Function,
  trian::GridapDistributed.DistributedTriangulation)

  fields = map(trian.trians) do t
    panelwise_cellfield(f,t)
  end
  GridapDistributed.DistributedCellField(fields,trian)
end

function geo_map_func(trian::DistributedTriangulation)
  println("distributed geo map")

  model = get_background_model(trian)
  owned_panel_ids = get_owned_panel_ids(model)
  return geo_map_func(owned_panel_ids)
end


function geo_map_func(owned_panel_ids::AbstractArray)
  println("distributed geo map")

  @assert typeof(owned_panel_ids) <: DebugArray || typeof(owned_panel_ids) <: MPIArray "\n Not distributed panel ids"

  cell_geo_map = map(owned_panel_ids) do pid
    return lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), pid)
  end
  return cell_geo_map
end
