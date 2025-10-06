function get_panel_ids(trian::GridapDistributed.DistributedTriangulation)
  gids = get_cell_gids(trian.model)

  dpanel_ids = map(trian.trians,partition(gids)) do trian, cid
    bgmodel = get_background_model(trian)
    panel_model = get_parent_model(bgmodel)
    panel_ids = get_panel_ids(panel_model)

    # return just the owned panel ids
    own2local = own_to_local(cid)
    local2global = local_to_global(cid)
    ids = local2global[own2local]

    return panel_ids[ids]
  end
  return dpanel_ids
end
