function get_panel_ids(trian::GridapDistributed.DistributedTriangulation)
  dpanel_ids = get_panel_ids(trian.model)
  panel_ids = map(trian.trians,dpanel_ids) do trian, panel_ids
    panel_ids = get_panel_ids(trian,panel_ids)
    return panel_ids
  end
  return panel_ids
end

function get_owned_panel_ids(trian::GridapDistributed.DistributedTriangulation)
  gids = get_cell_gids(trian.model)
  dpanel_ids = get_panel_ids(trian)

  owned_panel_ids = map(dpanel_ids,partition(gids)) do panel_ids, cids
    owned_cells = own_to_local(cids)
    return panel_ids[owned_cells]
  end
  return owned_panel_ids
end


function get_skel_panel_ids(skel_panel_ids::AbstractArray{<:SkeletonPair})
  skel_panel_ids_plus = map(skel_panel_ids) do s
    return s.plus
  end
  return skel_panel_ids_plus
end
