struct SphereAmbientCellMap{A,B,C,D} <: Map # map from ref FE -> panel p on sphere
  Rp1::A
  R1p::A
  Bump::B
  Gnomonic::C
  Sigma::D
end

function Gridap.Arrays.return_cache(f::SphereAmbientCellMap,panel_id::Int64,cmap)
  y = first(f.Rp1)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SphereAmbientCellMap,panel_id::Int64,cmap)
  y = cache
  y = f.R1p[panel_id] ∘ f.Sigma ∘ f.Gnomonic ∘ f.Bump ∘ f.Rp1[panel_id] ∘ cmap
  return y
end
