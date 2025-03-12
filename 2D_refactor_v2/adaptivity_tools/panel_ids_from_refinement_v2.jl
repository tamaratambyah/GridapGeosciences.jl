using Gridap

function get_panel_ids(model::DiscreteModel)
    return collect(1:num_cells(model))
end

function get_panel_ids(model::Gridap.Adaptivity.AdaptedDiscreteModel{Dc}) where Dc
  # println("no bang")
     panel_ids = copy(model.glue.n2o_faces_map[Dc+1])
     get_panel_ids!(panel_ids, model.parent)
     return panel_ids
end


function get_panel_ids!(panel_ids, model::Gridap.Adaptivity.AdaptedDiscreteModel{Dc}) where Dc
  # println("bang")
    n2o = model.glue.n2o_faces_map[Dc+1]
    panel_ids .= n2o[panel_ids]
    get_panel_ids!(panel_ids, model.parent)
end

function get_panel_ids!(panel_ids, model::DiscreteModel)
  # println("coarest model")
end


###############################################################################

function _get_panel_ids(model::Gridap.Adaptivity.AdaptedDiscreteModel{Dc}) where Dc

  panel_ids = copy(model.glue.n2o_faces_map[Dc+1])

  while Gridap.Adaptivity.is_child(model,model.parent) && (typeof(model.parent) <: Gridap.Adaptivity.AdaptedDiscreteModel)
    # println("do a thing")

    n2o = model.glue.n2o_faces_map[Dc+1]
    panel_ids .= n2o[panel_ids]
    # copy!(panel_ids,n2o[panel_ids])

    model = model.parent

  end

  return panel_ids
end
