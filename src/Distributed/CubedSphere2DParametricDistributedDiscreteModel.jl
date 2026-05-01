################################################################################
#### CubedSphereParametricDistributedDiscreteModel
#### Has serial CubedSphereParametricDiscreteModel underneath
#### The provided panel_ids are owned+ghost. Return only owned+ghost in get_panel_ids
#### Return owned in get_owned_panel_ids
#### Implemenet the interface of GridapDistributed.GenericDistributedDiscreteModel
################################################################################
# const CubedSphereParametricDistributedDiscreteModel{Dc,Dp,T} = GridapDistributed.GenericDistributedDiscreteModel{Dc,Dp,<:AbstractArray{T}} where T<:Union{<:CubedSphereParametricDiscreteModel{Dc,Dp},<:AdaptedDiscreteModel{Dc,Dp}}

const CubedSphere2DParametricDistributedDiscreteModel{T} = GridapDistributed.GenericDistributedDiscreteModel{2,2,<:AbstractArray{T}} where T<:Union{<:CubedSphere2DParametricDiscreteModel,<:AdaptedDiscreteModel{2,2}}
const CubedSphere3DParametricDistributedDiscreteModel{T} = GridapDistributed.GenericDistributedDiscreteModel{3,3,<:AbstractArray{T}} where T<:Union{<:CubedSphere3DParametricDiscreteModel,<:AdaptedDiscreteModel{3,3}}

const CubedSphereParametricDistributedDiscreteModel{T} = Union{CubedSphere2DParametricDistributedDiscreteModel{T},CubedSphere3DParametricDistributedDiscreteModel{T}}

# Note the dmodel passed in does not have CubedSphereParametricDiscreteModel as local models
# Thus, we need to explicitly pass the radius to the constructor as well
# I am unsure where to put the radius, store in metadata for now
function CubedSphere2DParametricDistributedDiscreteModel(
  model::GridapDistributed.DistributedDiscreteModel,
  dpanel_ids::AbstractArray,
  radius::Real
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
    # owned_cells = own_to_local(cids)
    # panel_ids = pids[owned_cells]
    CubedSphere2DParametricDiscreteModel(grid,topo,labels,pids,radius)
   end
  return GridapDistributed.GenericDistributedDiscreteModel(models,gids)
end

## these are local+ghost panel ids
## for dmodels where the local model is CubedSphereParametricDiscreteModel
function get_panel_ids(dmodel::CubedSphereParametricDistributedDiscreteModel{T}) where {T}
  return map(get_panel_ids,local_views(dmodel))
end

## these are local+ghost panel ids
## for omodels where the local model is AdaptedDiscreteModel
# function get_panel_ids(dmodel::CubedSphereParametricDistributedDiscreteModel{Dc,Dp,<:AdaptedDiscreteModel{Dc,Dp}}) where {Dc,Dp}
#   panel_ids = map(local_views(dmodel)) do lmodel
#     get_panel_ids(lmodel.model)
#   end
#   return panel_ids
# end

# return owned here to assist with triangulations
function get_owned_panel_ids(dmodel::CubedSphereParametricDistributedDiscreteModel{T}) where {T}
  gids = get_cell_gids(dmodel)
  dpanel_ids = get_panel_ids(dmodel)

  panel_ids = map(dpanel_ids,partition(gids)) do pids, cids
    owned_cells = own_to_local(cids)
    return pids[owned_cells]
  end
  return panel_ids

end

function get_forward_map_generator(dmodel::CubedSphereParametricDistributedDiscreteModel{T}) where {T}
  return map(get_forward_map_generator,local_views(dmodel))
end

function get_radius(dmodel::CubedSphereParametricDistributedDiscreteModel{T}) where {T}
  radii =  map(get_radius,local_views(dmodel))
  radius = zero(eltype(radii))
  map(radii) do r
    radius = r
  end
  return radius
end

function get_thickness(dmodel::CubedSphere2DParametricDistributedDiscreteModel{T}) where {T}
  @notimplemented """\n
  The model is two dimensional, get_thickness not defined.
  """
end

function get_thickness(dmodel::CubedSphere3DParametricDistributedDiscreteModel) where {T}
  Ts =  map(get_thickness,local_views(dmodel))
  thickness = zero(eltype(Ts))
  map(Ts) do t
    thickness = t
  end
  return thickness
end
