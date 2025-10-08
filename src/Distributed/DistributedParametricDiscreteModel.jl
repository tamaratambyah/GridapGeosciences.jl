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
    _grid = get_grid(model)

    ## construct the new grid by hand, so that specific cmaps are given
    nodes = collect1d(get_node_coordinates(_grid)) # these are just junk nodes, never used
    ctype = collect1d(get_cell_type(_grid))
    cmaps = collect1d(get_cell_map(_grid))
    grid = Gridap.Geometry.UnstructuredGrid(nodes,get_cell_node_ids(_grid),
              get_reffes(_grid),ctype,OrientationStyle(_grid),
                        nothing,cmaps)


    topo = UnstructuredGridTopology(get_grid_topology(model))
    labels = FaceLabeling(topo)

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
