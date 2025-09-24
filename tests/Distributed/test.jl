using DrWatson
using Gridap
using GridapDistributed
using PartitionedArrays
using MPIPreferences

using Gridap.Geometry
using Gridap.Adaptivity
MPIPreferences.use_jll_binary()

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
s_model_coarse = coarse_parametric_model()
s_model_ref = refine(s_model_coarse)
s_model_ref_ref = refine(s_model_ref)

spanel_ids = [get_panel_ids(s_model_ref_ref),get_panel_ids(s_model_ref),get_panel_ids(s_model_coarse)]

# s_model_coarse = Geometry.UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(2,3)))
# s_model_ref = refine(s_model_coarse)
# s_model_ref_ref = refine(s_model_ref)

################################################################################
##### Distributed models
################################################################################
nprocs = 4
ranks  = with_debug() do distribute
  distribute(LinearIndices((4,)))
end


part_to_cells = [PartitionedArrays.local_range(rank,nprocs,num_cells(s_model_coarse)) for rank in 1:nprocs]
coarse_cell_to_part = zeros(Int32,num_cells(s_model_coarse))
for (rank, cells) in enumerate(part_to_cells)
  coarse_cell_to_part[cells] .= rank
end



models = [s_model_ref_ref.model,s_model_ref.model,s_model_coarse]
glues = [s_model_ref_ref.glue,s_model_ref.glue]
cell_to_part = Vector{Any}(undef,length(models))
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



dmodels = Vector{Any}(undef,length(models))
dpanel_ids = Vector{Any}(undef,length(models))
owned_panel_ids = Vector{Any}(undef,length(models))

dmodels[end] = DiscreteModel(ranks,models[end],cell_to_part[end])
dpanel_ids[end],owned_panel_ids[end] = distributed_panel_ids(dmodels[end],spanel_ids[end])
for level in length(models)-1:-1:1
  child = Gridap.Geometry.UnstructuredDiscreteModel(DiscreteModel(ranks,models[level],cell_to_part[level]))
  parent = dmodels[level+1]
  glue = DistributedAdaptivityGlue(glues[level],parent,child)
  dmodels[level] = GridapDistributed.DistributedAdaptedDiscreteModel(child,parent,glue)
  dpanel_ids[level],owned_panel_ids[level] = distributed_panel_ids(child,spanel_ids[level])
end

for i in 1:length(models)
  println("level ", i)
  map(local_views(dpanel_ids[i]),local_views(owned_panel_ids[i])) do p,op
    println("Panel ids: ", p)
    println("Owned ids: ", op)
  end
end


include("vtk.jl")
for (level,(_dmodel,own_panel_ids)) in enumerate(zip(dmodels,owned_panel_ids))
  dmodel = (level == 3) ? _dmodel : Gridap.Adaptivity.get_model(_dmodel)

  # gids = get_cell_gids(dmodel)
  cell_geo_map = map(own_panel_ids) do panel_ids
    return lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
  end

  writevtk(Triangulation(dmodel),dir*"/ambient_model_ref_$(level)"; append=false)#, geo_map=cell_geo_map)
end
