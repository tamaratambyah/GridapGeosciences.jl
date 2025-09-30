using DrWatson
using Gridap
using GridapDistributed
using PartitionedArrays
using MPIPreferences

using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers
MPIPreferences.use_jll_binary()

dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


### This is Jordi's function from GridapDistributed/test/AdaptivityTests.jl
function DistributedAdaptivityGlue(serial_glue,parent,child)
  glue = map(partition(get_cell_gids(parent)),partition(get_cell_gids(child))) do parent_gids, child_gids
    old_g2l = global_to_local(parent_gids)
    old_l2g = local_to_global(parent_gids)
    new_l2g = local_to_global(child_gids)

    n2o_cell_map  = lazy_map(Reindex(old_g2l),serial_glue.n2o_faces_map[3][new_l2g])
    n2o_faces_map = [Int64[],Int64[],collect(n2o_cell_map)]
    n2o_cell_to_child_id = serial_glue.n2o_cell_to_child_id[new_l2g]
    rrules = serial_glue.refinement_rules[old_l2g]
    Gridap.Adaptivity.AdaptivityGlue(n2o_faces_map,n2o_cell_to_child_id,rrules)
  end
  return glue
end


#### The distributed panel ids are extracted from the serial. This includes both
#### owned+ghost panel_ids.
#### The owned panel_ids are extracted by determining the owned cell
function distributed_panel_ids(dmodel,spanel_ids::AbstractArray{Int})
  gids = get_cell_gids(dmodel)

  dpanel_ids = map(partition(gids)) do ids
    lid_to_gid = local_to_global(ids)
    return spanel_ids[lid_to_gid]
  end

  owned_panel_ids = map(dpanel_ids,partition(gids)) do panel_ids, cids
    owned_cells = own_to_local(cids)
    return panel_ids[owned_cells]
  end

  return dpanel_ids, owned_panel_ids
end

################################################################################
using GridapGeosciences
n_ref_lvls = 2
s_models  = get_refined_models(n_ref_lvls,true)
spanel_ids = map(m->get_panel_ids(m),s_models)
s_model_coarse = s_models[end]



################################################################################
##### Distributed models
################################################################################
nprocs = 6
ranks  = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end


part_to_cells = [PartitionedArrays.local_range(rank,nprocs,num_cells(s_model_coarse)) for rank in 1:nprocs]
coarse_cell_to_part = zeros(Int32,num_cells(s_model_coarse))
for (rank, cells) in enumerate(part_to_cells)
  coarse_cell_to_part[cells] .= rank
end


models = map(m->get_model(m),s_models[1:end-1])
glues = map(m->get_adaptivity_glue(m),s_models[1:end-1])
push!(models,s_model_coarse)


cell_to_part = Vector{Vector{Int32}}(undef,length(models))
cell_to_part[end] = coarse_cell_to_part
for level in length(models)-1:-1:1
  n2o_cells = glues[level].n2o_faces_map[3]
  cell_to_part[level] = cell_to_part[level+1][n2o_cells]
end

## plot the partition in serial
for (level,(model,cparts,panel_ids)) in enumerate(zip(models,cell_to_part,spanel_ids))
  geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
  writevtk(Triangulation(model),dir*"/serial_model_ref_$(level)";append=false,celldata=["part" => cparts],geo_map=geo_map)
end

dmodels = Vector{GridapDistributed.DistributedDiscreteModel}(undef,length(models))
dpanel_ids = Vector{AbstractArray{Vector{Int}}}(undef,length(models))
owned_panel_ids = Vector{AbstractArray{Vector{Int}}}(undef,length(models))

dmodels[end] = DiscreteModel(ranks,models[end],cell_to_part[end])
dpanel_ids[end],owned_panel_ids[end] = distributed_panel_ids(dmodels[end],spanel_ids[end])
for level in length(models)-1:-1:1
  child = DiscreteModel(ranks,models[level],cell_to_part[level])
  parent = dmodels[level+1]
  glue = DistributedAdaptivityGlue(glues[level],parent,child)
  dmodels[level] = GridapDistributed.DistributedAdaptedDiscreteModel(child,parent,glue)
  dpanel_ids[level],owned_panel_ids[level] = distributed_panel_ids(child,spanel_ids[level])
end

for i in 1:length(models)
  println("Ref lvl: ", length(models)-i)
  map(local_views(dpanel_ids[i]),local_views(owned_panel_ids[i]),ranks) do p,op, r
    println("Proc no: $r")
    println("\tOwned+Ghost panel_ids: ", p)
    println("\tOwned panel_ids: ", op)
  end
end

using GridapSolvers

mh = GridapSolvers.ModelHierarchy(dmodels)

typeof(dmodels) <:Vector{<:GridapDistributed.DistributedDiscreteModel}

function ModelHierarchy(models::Vector{<:GridapDistributed.DistributedDiscreteModel})
  nlevs = length(models)
  # @check all(map(i -> isa(models[i],GridapDistributed.DistributedAdaptedDiscreteModel),1:nlevs-1)) "Hierarchy models are not adapted models."
  @check all(map(i -> isa(models[i],DistributedAdaptedParametricDiscreteModel),1:nlevs-1)) "Hierarchy models are not adapted models."
  for lev in 1:nlevs-1
    @check Adaptivity.is_child(models[lev],models[lev+1]) "Incorrect hierarchy of models."
  end
  ranks = get_parts(models[1])
  @check all(m -> length(get_parts(m)) === length(ranks), models) "Models have different communicators."

  level_parts = fill(ranks,nlevs)
  meshes = Vector{GridapSolvers.MultilevelTools.ModelHierarchyLevel}(undef,nlevs)
  for lev in 1:nlevs-1
    glue  = Gridap.Adaptivity.get_adaptivity_glue(models[lev])
    model = models[lev]
    meshes[lev] = GridapSolvers.MultilevelTools.ModelHierarchyLevel(lev,model,glue,nothing,nothing)
  end
  meshes[nlevs] = GridapSolvers.MultilevelTools.ModelHierarchyLevel(nlevs,models[nlevs],nothing,nothing,nothing)
  return GridapSolvers.MultilevelTools.HierarchicalArray(meshes,level_parts)
end


1;



#### Plot models 1 at a time. If using a loop, vtk crashes
include("vtk.jl")

level = 1
dmodel = dmodels[level]
o_pids = owned_panel_ids[level]

cell_geo_map = map(o_pids) do pid
  return lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), pid)
end

writevtk(Triangulation(dmodel),dir*"/ambient_model_ref_$(level)", append=false, compress=false,geo_map=cell_geo_map)



map(local_views(dmodel),ranks,o_pids) do model, r, pid
  println(typeof(model.model))
  grid = get_grid(model)

  cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), pid)
  writevtk(Triangulation(model.model),dir*"/ambient_model_ref_1_rank$r", append=false, compress=false,geo_map=cell_geo_map)

end
