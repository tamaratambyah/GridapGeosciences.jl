function get_panel_ids(trian::GridapDistributed.DistributedTriangulation{Df}) where Df

  panel_ids = map(trian.trians) do t
    return get_panel_ids(t)
  end
  return panel_ids
end

function get_skel_panel_ids(skel_panel_ids::AbstractArray{<:SkeletonPair})
  skel_panel_ids_plus = map(skel_panel_ids) do s
    return s.plus
  end
  return skel_panel_ids_plus
end
