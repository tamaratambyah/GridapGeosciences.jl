function get_panel_ids(strian::SkeletonTriangulation)
  plus = get_face_panel_ids(strian.plus)
  minus = get_face_panel_ids(strian.plus)
  SkeletonPair(plus,minus)
end

function Geometry.get_facet_normal(trian::SkeletonTriangulation,cell_geo_map::AbstractArray)
  println("skeleton facet normal")
  plus = get_facet_normal(trian.plus,cell_geo_map)
  minus = get_facet_normal(trian.minus,cell_geo_map)
  SkeletonPair(plus,minus)
end

function pushforward_normal(trian::SkeletonTriangulation)
  plus, J_plus = pushforward_normal(trian.plus)
  minus, J_minus = pushforward_normal(trian.minus)
  SkeletonPair(plus,minus), SkeletonPair(J_plus,J_minus)
end


function pullback_area_form(trian::SkeletonTriangulation)
  plus = pullback_area_form(trian.plus)
  minus = pullback_area_form(trian.minus)
  SkeletonPair(plus,minus)
end
