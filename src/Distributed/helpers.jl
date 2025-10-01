

#### The distributed panel ids are extracted from the serial. This includes both
#### owned+ghost panel_ids.
#### The owned panel_ids are extracted by determining the owned cell
function distributed_panel_ids(dmodel,spanel_ids::AbstractArray{Int})
  gids = get_cell_gids(dmodel)

  dpanel_ids = map(partition(gids)) do ids
    lid_to_gid = local_to_global(ids)
    return spanel_ids[lid_to_gid]
  end

  owned_panel_ids = map(dpanel_ids,partition(gids)) do panel_ids, cids
    owned_cells = own_to_local(cids)
    return panel_ids[owned_cells]
  end

  return dpanel_ids, owned_panel_ids
end
