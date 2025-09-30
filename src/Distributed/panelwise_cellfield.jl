function panelwise_cellfield(f::Function,
  trian::GridapDistributed.DistributedTriangulation,
  panel_ids::AbstractArray)

  @assert typeof(panel_ids) <: DebugArray || typeof(panel_ids) <: MPIArray "\n Not distributed panel ids"

  fields = map(trian.trians,panel_ids) do trian, pids
    panelwise_cellfield(f,trian,pids)
  end
  GridapDistributed.DistributedCellField(fields,trian)
end


function geo_map_func(owned_panel_ids::AbstractArray)
  println("distributed geo map")

  @assert typeof(owned_panel_ids) <: DebugArray || typeof(owned_panel_ids) <: MPIArray "\n Not distributed panel ids"

  cell_geo_map = map(owned_panel_ids) do pid
    return lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), pid)
  end
  return cell_geo_map
end
