"""
get_panel_ids

returns the panel id = 1,…,6, for unrefined and refined cubed models
It is assumed that the coarset model has 1 cell per panel
"""
function get_panel_ids(args...)
  @abstractmethod
end

function get_panel_ids(model::DiscreteModel)
  return collect(1:num_cells(model))
end

function get_panel_ids!(panel_ids, model::DiscreteModel)
end


function get_panel_ids(model::Gridap.Adaptivity.AdaptedDiscreteModel{Dc}) where Dc
  panel_ids = copy(model.glue.n2o_faces_map[Dc+1])
  get_panel_ids!(panel_ids, model.parent)
  return panel_ids
end


function get_panel_ids!(panel_ids, model::Gridap.Adaptivity.AdaptedDiscreteModel{Dc}) where Dc
  n2o = model.glue.n2o_faces_map[Dc+1]
  panel_ids .= n2o[panel_ids]
  get_panel_ids!(panel_ids, model.parent)
end


function geo_map_func(trian::Triangulation)
  panel_ids = get_panel_ids(trian)
  geo_map_func(panel_ids)
end

function geo_map_func(panel_ids::AbstractArray{Int})
  # println("serial geo map")
  return lazy_map(p -> ForwardMap(p), panel_ids)
end

### latlon geo func
function latlon_geo_map_func(trian::Triangulation)
  panel_ids = get_panel_ids(trian)
  latlon_geo_map_func(panel_ids)
end

### here we have to compose separate maps so vtk uses the cellwise-version of
### Cartesian2SphereicalMap()
function latlon_geo_map_func(panel_ids::AbstractArray{Int})
  # println("latolon serial geo map")

  cell_geo_map = geo_map_func(panel_ids)
  fi = lazy_map(p->Cartesian2SphereicalMap(),panel_ids)
  latlon_cell_geo_map = lazy_map(∘, fi, cell_geo_map)

  return latlon_cell_geo_map
end
