function pullback_area_form(trian::TriangulationView)
  println("trian view")
  _pullback_area_form(trian.parent)
end

function Gridap.Geometry.get_facet_normal(trian::Gridap.Geometry.TriangulationView,cell_geo_map::AbstractArray)
  Gridap.Geometry.get_facet_normal(trian.parent,cell_geo_map)
end
