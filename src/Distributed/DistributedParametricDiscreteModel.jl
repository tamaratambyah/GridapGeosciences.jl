################################################################################
#### DistributedParametricDiscreteModel
#### The panel_ids are owned+ghost. Return only the owned in get_panel_ids
#### Implemenet the interface of GridapDistributed.GenericDistributedDiscreteModel
################################################################################
struct DistributedParametricDiscreteModel{Dc,Dp,A,B,C,D} <: GridapDistributed.DistributedDiscreteModel{Dc,Dp}
  models::A
  face_gids::B
  panel_ids::C
  metadata::D
end

function DistributedParametricDiscreteModel(models::AbstractArray{<:DiscreteModel{Dc,Dp}},
  gids::PRange, panel_ids::AbstractArray;
  metadata=nothing) where {Dc,Dp}

  face_gids=Vector{PRange}(undef,Dc+1)
  face_gids[Dc+1] = gids

  A = typeof(models)
  B = typeof(face_gids)
  C = typeof(panel_ids)
  D = typeof(metadata)
  DistributedParametricDiscreteModel{Dc,Dp,A,B,C,D}(models,face_gids,panel_ids,metadata)
end

function DistributedParametricDiscreteModel(dmodel::GridapDistributed.GenericDistributedDiscreteModel{Dc,Dp},panel_ids) where {Dc,Dp}
  models = dmodel.models
  face_gids = dmodel.face_gids
  metadata = dmodel.metadata
  A = typeof(models)
  B = typeof(face_gids)
  C = typeof(panel_ids)
  D = typeof(metadata)

  DistributedParametricDiscreteModel{Dc,Dp,A,B,C,D}(models,face_gids,panel_ids,metadata)
end

function get_panel_ids(dmodel::DistributedParametricDiscreteModel)
  dpanel_ids = dmodel.panel_ids
  gids = get_cell_gids(dmodel)

  owned_panel_ids = map(dpanel_ids,partition(gids)) do panel_ids, cids
    owned_cells = own_to_local(cids)
    return panel_ids[owned_cells]
  end

  return owned_panel_ids
end

### The remaining interface is the same as GridapDistributed.GenericDistributedDiscreteModel
GridapDistributed.local_views(a::DistributedParametricDiscreteModel) = a.models

function GridapDistributed.get_cell_gids(model::DistributedParametricDiscreteModel{Dc}) where Dc
  model.face_gids[Dc+1]
end

function GridapDistributed.get_face_gids(model::DistributedParametricDiscreteModel,dim::Integer)
  GridapDistributed._setup_face_gids!(model,dim)
  model.face_gids[dim+1]
end

function GridapDistributed._setup_face_gids!(dmodel::DistributedParametricDiscreteModel{Dc},dim) where {Dc}
  Gridap.Helpers.@check 0 <= dim <= Dc
  if !isassigned(dmodel.face_gids,dim+1)
    mgids   = dmodel.face_gids[Dc+1]
    nlfaces = map(local_views(dmodel)) do model
      num_faces(model,dim)
    end
    cell_lfaces = map(local_views(dmodel)) do model
      topo  = get_grid_topology(model)
      faces = get_faces(topo, Dc, dim)
    end
    dmodel.face_gids[dim+1] = generate_gids(mgids,cell_lfaces,nlfaces)
  end
  return
end




################################################################################
#### DistributedAdaptedParametricDiscreteModel
#### The panel_ids are owned+ghost
#### Implemenet the interface of GridapDistributed.DistributedAdaptedDiscreteModel
################################################################################
const DistributedAdaptedParametricDiscreteModel{Dc,Dp} = DistributedParametricDiscreteModel{Dc,Dp,<:AbstractArray{<:AdaptedDiscreteModel{Dc,Dp}}}

function DistributedAdaptedParametricDiscreteModel(
  model  :: GridapDistributed.DistributedDiscreteModel,
  parent :: GridapDistributed.DistributedDiscreteModel,
  glue   :: AbstractArray{<:AdaptivityGlue},
  panel_ids::AbstractArray
)
  models = map(local_views(model),local_views(parent),glue) do model, parent, glue
    AdaptedDiscreteModel(model,parent,glue)
  end
  gids = get_cell_gids(model)
  metadata = hasproperty(model,:metadata) ? model.metadata : nothing
  return DistributedParametricDiscreteModel(models,gids,panel_ids;metadata)
end

function Adaptivity.get_model(model::DistributedAdaptedParametricDiscreteModel)
  DistributedParametricDiscreteModel(
    map(get_model,local_views(model)),
    get_cell_gids(model),
    get_panel_ids(model);
    metadata = hasproperty(model,:metadata) ? model.metadata : nothing
  )
end

function Adaptivity.get_parent(model::DistributedAdaptedParametricDiscreteModel)
  msg = " Error: Cannot get global parent model. \n
          We do not keep the global ids of the parent model within the children.\n
          You can extract the local parents with map(get_parent,local_views(model))"
  @notimplemented msg
end

function Adaptivity.get_adaptivity_glue(model::DistributedAdaptedParametricDiscreteModel)
  return map(Adaptivity.get_adaptivity_glue,local_views(model))
end
