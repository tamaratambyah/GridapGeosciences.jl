function get_panel_ids(trian::GridapDistributed.DistributedTriangulation)
  # dpanel_ids = map(trian.trians) do trian
  #   panel_ids = get_panel_ids(trian)
  #   return panel_ids
  # end
  dpanel_ids = trian.model.metadata.dpanel_ids

  dpanel_ids = map(trian.trians,dpanel_ids) do trian, panel_ids
    panel_ids = get_panel_ids(trian,panel_ids)
    return panel_ids
  end
  return dpanel_ids
end
