############ debug
models = get_refined_models(2)
panel_model = models[1]
p_fe = 1

panel_ids = get_panel_ids(panel_model)
degree = 2*(p_fe + 1)

Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,degree)

u_cf = panelwise_cellfield(u,Ω_panel,panel_ids)

using Gridap.Helpers, Gridap.Geometry, Gridap.CellData, Gridap.Fields
function GridapGeosciences.panelwise_cellfield(f::Function,trian::BodyFittedTriangulation,panel_ids::AbstractArray{Int})
  println("new cell fields")
  @check length(panel_ids) == num_cells(trian) "\n Incorrect panel ids"
  cell_field = map(p->GenericField(f(p)),panel_ids)
  CellData.GenericCellField(cell_field,trian,PhysicalDomain())

  _cf = _cell_data(f,trian,panel_ids)
  CellData.GenericCellField(_cf,trian,ReferenceDomain())

end

function _cell_data(f::Function,trian,panel_ids::AbstractArray)
  # make physical cf
  cell_field = map(p->GenericField(f(p)),panel_ids)
  cf = CellData.GenericCellField(cell_field,trian,PhysicalDomain())

  panel_model = get_background_model(trian)
  cmap = get_cell_map(get_grid(panel_model))

  # compose with the map from refernce cell -> physical cell
  lazy_map(∘,get_data(cf),cmap)

end
