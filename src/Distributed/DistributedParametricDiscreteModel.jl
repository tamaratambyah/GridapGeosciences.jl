################################################################################
#### DistributedParametricDiscreteModel
#### Has serial ParametricDiscreteModel underneath
#### The provided panel_ids are owned+ghost. Return only the owned in get_panel_ids
#### Implemenet the interface of GridapDistributed.GenericDistributedDiscreteModel
################################################################################
const DistributedParametricDiscreteModel{Dc,Dp} = GridapDistributed.GenericDistributedDiscreteModel{Dc,Dp,<:AbstractArray{<:ParametricDiscreteModel{Dc,Dp}}}

struct DistributedParametricDiscreteModelCache{A}
  dpanel_ids::A
end

function DistributedParametricDiscreteModel(
  model::GridapDistributed.DistributedDiscreteModel,
  dpanel_ids::AbstractArray
)
  gids = get_cell_gids(model)

  models = map(local_views(model),dpanel_ids,partition(gids)) do model, pids, cids
    grid = get_grid(model)
    topo = get_grid_topology(model)
    labels = get_face_labeling(model)

    # extract the owned panel ids
    owned_cells = own_to_local(cids)
    panel_ids = pids[owned_cells]
    ParametricDiscreteModel(grid,topo,labels,panel_ids)
   end

  metadata = DistributedParametricDiscreteModelCache(dpanel_ids)
  return GridapDistributed.GenericDistributedDiscreteModel(models,gids;metadata)
end

function get_panel_ids(dmodel::DistributedParametricDiscreteModel)
  return map(get_panel_ids,local_views(dmodel))
end
