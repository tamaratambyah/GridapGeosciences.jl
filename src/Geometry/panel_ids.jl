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
