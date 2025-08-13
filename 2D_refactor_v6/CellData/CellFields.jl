function panelwise_cellfield(f::Function,trian::Triangulation,panel_ids::AbstractArray{Int})
  cell_field = map(p->GenericField(f(p)),panel_ids)
  CellData.GenericCellField(cell_field,trian,PhysicalDomain())
end

function ambient_cellfield(panel_cf::CellField,ambient_trian::Triangulation,panel_ids::AbstractArray{Int})
  inv_f = lazy_map(p->InverseMap(p),panel_ids)
  _cf = change_domain(panel_cf,DomainStyle(panel_cf),PhysicalDomain())
  cf_mapped = lazy_map(Broadcasting(∘),get_data(_cf),inv_f)
  CellData.GenericCellField(cf_mapped,ambient_trian,PhysicalDomain() )
end
