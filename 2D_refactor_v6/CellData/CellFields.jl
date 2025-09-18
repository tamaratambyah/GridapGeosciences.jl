function panelwise_cellfield(f::Function,trian::Triangulation,panel_ids::AbstractArray{Int})
  @check length(panel_ids) == num_cells(trian) "\n Incorrect panel ids"
  cell_field = map(p->GenericField(f(p)),panel_ids)
  CellData.GenericCellField(cell_field,trian,PhysicalDomain())
end

### For Adapted Triangulations, return cellfield on the atrian
function panelwise_cellfield(f::Function,atrian::AdaptedTriangulation)
  println("Adapted Boundary Triangulation")
  cf = panelwise_cellfield(f,atrian.trian)
  panelwise_cellfield(cf,atrian)
end

function panelwise_cellfield(cf::CellField,atrian::AdaptedTriangulation)
  GenericCellField(get_data(cf),atrian,DomainStyle(cf))
end

function panelwise_cellfield(cf::SkeletonPair,atrian::AdaptedTriangulation)
  plus = GenericCellField(get_data(cf.plus),atrian,ReferenceDomain())
  minus = GenericCellField(get_data(cf.minus),atrian,ReferenceDomain())
  return SkeletonPair(plus,minus)
end

### For Skeleton and Boundary Triangulations, map with map from reference face
### Thus, return cellfield on the reference domain

### Boundary
function panelwise_cellfield(f::Function,trian::BoundaryTriangulation)
  _face_cf = boundary_cell_data(f,trian)
  face_cf = GenericCellField(_face_cf,trian,ReferenceDomain())
  return face_cf
end

function boundary_cell_data(f::Function,trian::BoundaryTriangulation)
  face_panel_ids = get_panel_ids(trian)

  @check length(face_panel_ids) == num_cells(trian) "\n Incorrect panel ids"

  # make physical cf
  cell_field = map(p->GenericField(f(p)),face_panel_ids)
  cf = CellData.GenericCellField(cell_field,trian,PhysicalDomain())

  glue = get_glue(trian,Val(2))

  panel_model = get_background_model(trian)
  cmap = get_cell_map(get_grid(panel_model))
  fcmap = lazy_map(Reindex(cmap),glue.tface_to_mface)

  # compose with the map from reference face -> refernce cell -> physical cell
  ref_face_2_phys_cell_map = lazy_map(∘,fcmap,glue.tface_to_mface_map)
  face_cf = lazy_map(∘,get_data(cf),ref_face_2_phys_cell_map)
  return face_cf
end

### Skeleton - left and right boundary triangulations returned as SkeletonPair
function panelwise_cellfield(f::Function,trian::SkeletonTriangulation)
  ### plus
  _face_cf_plus = boundary_cell_data(f,trian.plus)
  plus = GenericCellField(_face_cf_plus,trian,ReferenceDomain())

  #### minus
  _face_cf_minus = boundary_cell_data(f,trian.minus)
  minus = GenericCellField(_face_cf_minus,trian,ReferenceDomain())

  SkeletonPair(plus,minus)
end






function ambient_cellfield(panel_cf::CellField,ambient_trian::Triangulation,panel_ids::AbstractArray{Int})
  inv_f = lazy_map(p->InverseMap(p),panel_ids)
  _cf = change_domain(panel_cf,DomainStyle(panel_cf),PhysicalDomain())
  cf_mapped = lazy_map(Broadcasting(∘),get_data(_cf),inv_f)
  CellData.GenericCellField(cf_mapped,ambient_trian,PhysicalDomain() )
end
