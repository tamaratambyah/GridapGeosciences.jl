struct CubeParametricCellMap{A,B} <: Map # map from ref FE -> panel p of cube
  Rp1::A
  R1p::A
  Bump::B
end

function Gridap.Arrays.return_cache(f::CubeParametricCellMap,panel_id::Int64,cmap)
  y = first(f.Rp1)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::CubeParametricCellMap,panel_id::Int64,cmap)
  y = cache
  # y = f.Bump ∘ f.Rp1[panel_id] ∘ cmap
  y = Operation(f.Bump)(Operation(f.Rp1[panel_id])(cmap) )

  # About the same speed. ∘ is easier to read
  return y
end

_model = cube_model_3D
panel_ids = get_panel_ids(_model)
Ω = Triangulation(_model)

R(x) = rp1[x]
# _R = Fill(  CellField(R(panel_ids),Ω,PhysicalDomain()), 6)


_R = [CellField(rp1[i],Ω,PhysicalDomain()) for i in panel_ids]

cmap = get_cell_map(_model)

_cmap = lazy_map(∘,_R,cmap)

_cell_Jt = lazy_map(∇,_cmap)
_cell_Jtx = lazy_map(evaluate,_cell_Jt,quad.cell_point)
