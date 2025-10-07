function get_panel_ids(trian::GridapDistributed.DistributedTriangulation)
  dpanel_ids = map(trian.trians) do trian
    panel_ids = get_panel_ids(trian)
    return panel_ids
  end

  return dpanel_ids
end
