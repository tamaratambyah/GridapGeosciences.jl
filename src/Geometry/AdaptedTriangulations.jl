
function pushforward_normal(trian::AdaptedTriangulation,cell_geo_map::AbstractArray)
  cell_vectors = get_facet_normal(trian,cell_geo_map)
  get_normal_vector(trian,cell_vectors)
end

Geometry.get_facet_normal(trian::AdaptedTriangulation,cell_geo_map::AbstractArray) = Geometry.get_facet_normal(trian.trian,cell_geo_map)

function pushforward_normal(trian::AdaptedTriangulation)
  pushforward_normal(trian.trian)
end

function pullback_area_form(atrian::AdaptedTriangulation)
  cf = pullback_area_form(atrian.trian)
  ParametricCellField(cf,atrian)
end
