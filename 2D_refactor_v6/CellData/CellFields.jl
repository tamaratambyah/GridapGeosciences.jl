function panelwise_cellfield(f::Function,trian::Triangulation,panel_ids::AbstractArray{Int})
  cell_field = map(p->GenericField(f(p)),panel_ids)
  CellData.GenericCellField(cell_field,trian,PhysicalDomain())
end



function panelwise_cellfield(f::Function,atrian::AdaptedTriangulation,panel_ids::SkeletonPair{<:AbstractArray{Int}})
  panelwise_cellfield(f,atrian.trian,panel_ids)
end

function panelwise_cellfield(f::Function,trian::Triangulation,skel_panel_ids::SkeletonPair{<:AbstractArray{Int}})
  plus = panelwise_cellfield(f, trian.plus, skel_panel_ids.plus)
  minus = panelwise_cellfield(f, trian.minus, skel_panel_ids.minus)
  SkeletonPair(plus,minus)
end


function ambient_cellfield(panel_cf::CellField,ambient_trian::Triangulation,panel_ids::AbstractArray{Int})
  inv_f = lazy_map(p->InverseMap(p),panel_ids)
  _cf = change_domain(panel_cf,DomainStyle(panel_cf),PhysicalDomain())
  cf_mapped = lazy_map(Broadcasting(∘),get_data(_cf),inv_f)
  CellData.GenericCellField(cf_mapped,ambient_trian,PhysicalDomain() )
end
