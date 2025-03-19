struct CubePhysCellMap{A,B} <: Map # map from ref FE -> panel p of cube
  Rp1::A
  R1p::A
  Bump::B
end

function Gridap.Arrays.return_cache(f::CubePhysCellMap,panel_id::Int64,cmap)
  y = first(f.Rp1)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::CubePhysCellMap,panel_id::Int64,cmap)
  y = cache
  # y = f.R1p[panel_id] ∘ f.Bump∘ f.Bump∘ f.Rp1[panel_id] ∘ cmap
  y = f.Bump∘ f.Rp1[panel_id] ∘ cmap
  # y = Operation(f.R1p[panel_id])(Operation(f.Bump)(Operation(f.Bump)(Operation(f.Rp1[panel_id])(cmap) )  ) )

  # About the same speed. ∘ is easier to read
  return y
end
