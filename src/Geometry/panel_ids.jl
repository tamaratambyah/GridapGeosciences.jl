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

# need this global_pids junk array for dispatching in Distributed
function get_panel_ids(trian::BodyFittedTriangulation,global_pids=[1])
  panel_model = get_background_model(trian)
  @check typeof(panel_model) <: ParametricDiscreteModel "\n Not a ParametricDiscreteModel"
  get_panel_ids(panel_model)
end

function get_panel_ids(trian::Gridap.Geometry.TriangulationView,panel_ids::AbstractArray)
  println("new panel ids")
  panel_model = get_background_model(trian)
  Dc = num_cell_dims(panel_model)
  glue = get_glue(trian,Val(Dc))
  face_2_cell = glue.tface_to_mface
  face_panel_ids = panel_ids[face_2_cell]
  return face_panel_ids
end

function geo_map_func(panel_ids::AbstractArray{Int})
  println("serial geo map")
  return lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
end
