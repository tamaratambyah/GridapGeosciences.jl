function get_panel_ids(atrian::AdaptedTriangulation)
  get_panel_ids(atrian.trian)
end

Geometry.get_facet_normal(trian::AdaptedTriangulation,cell_geo_map::AbstractArray) = Geometry.get_facet_normal(trian.trian,cell_geo_map)

function pushforward_normal(trian::AdaptedTriangulation)
  pushforward_normal(trian.trian)
end
