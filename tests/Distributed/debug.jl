using Gridap
using GridapGeosciences
using DrWatson
using FillArrays
using Gridap.Geometry, Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Test

@testset "PanelIds" begin include("../GeometryTests/PanelIdsTests.jl") end



using GridapDistributed
using PartitionedArrays
using MPIPreferences
MPIPreferences.use_jll_binary()



dir = datadir("Distributed")
include("../convergence_tools.jl")


n_ref_lvls = 2

models  = get_refined_models(n_ref_lvls,true)


nprocs = 6

ranks = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end


s_models  = get_refined_models(n_ref_lvls,true)
spanel_ids = map(m->get_panel_ids(m),s_models)

s_model_coarse = s_models[end]

# extract the models and glues in arrays
models = map(m->Adaptivity.get_model(m),s_models[1:end-1])
glues = map(m->get_adaptivity_glue(m),s_models[1:end-1])
if typeof(s_model_coarse) <: AdaptedDiscreteModel
  push!(models,Adaptivity.get_model(s_model_coarse))
else
  push!(models,s_model_coarse)
end

# partition the processors
part_to_cells = [PartitionedArrays.local_range(rank,nprocs,num_cells(s_model_coarse)) for rank in 1:nprocs]

# store the partition for each model
cell_to_part = Vector{Vector{Int32}}(undef,length(models))

# get the coarse partition
coarse_cell_to_part = zeros(Int32,num_cells(s_model_coarse))
for (rank, cells) in enumerate(part_to_cells)
  coarse_cell_to_part[cells] .= rank
end
cell_to_part[end] = coarse_cell_to_part

# get the refine partition based on glue
for level in length(models)-1:-1:1
  n2o_cells = glues[level].n2o_faces_map[3]
  cell_to_part[level] = cell_to_part[level+1][n2o_cells]
end

# construct array of distributed models, distributed panel ids (all + owned only)
dmodels = Vector{DistributedParametricDiscreteModel}(undef,length(models))
# dmodels = Vector{Any}(undef,length(models))
dpanel_ids = Vector{AbstractArray{Vector{Int}}}(undef,length(models))
owned_panel_ids = Vector{AbstractArray{Vector{Int}}}(undef,length(models))

# the coarsest model is the last in the list
coarse_dmodel = DiscreteModel(ranks,models[end],cell_to_part[end])
dpanel_ids[end],owned_panel_ids[end] = distributed_panel_ids(coarse_dmodel,spanel_ids[end])
dmodels[end] = DistributedParametricDiscreteModel(coarse_dmodel,dpanel_ids[end])

# loop backwards through refinement levels
for level in length(models)-1:-1:1
  child = DiscreteModel(ranks,models[level],cell_to_part[level])
  parent = dmodels[level+1]
  glue = DistributedAdaptivityGlue(glues[level],parent,child)
  dpanel_ids[level],owned_panel_ids[level] = distributed_panel_ids(child,spanel_ids[level])

  dmodels[level] = DistributedParametricDiscreteModel(child, dpanel_ids[level])
end

################################################################################
#### BodyFittedTriangulation
################################################################################


dpanel_model = dmodels[2]
get_panel_ids(dpanel_model)

trian = Triangulation(dpanel_model)
get_panel_ids(trian)
get_owned_panel_ids(trian)
cell_geo_map = geo_map_func(get_owned_panel_ids(trian))
writevtk(trian,dir*"/distributed_model",append=false,geo_map=cell_geo_map)

################################################################################
#### BoundaryTriangulation
#### Need to return the mask that is all the interior cells
################################################################################

btrian = BoundaryTriangulation(dpanel_model)
b_panel_ids =  get_panel_ids(btrian)

b_geo_map = geo_map_func(b_panel_ids)
writevtk(btrian,dir*"/distributed_boundary_trian",append=false,geo_map=b_geo_map)



################################################################################
#### SkeletonTriangulation
#### Need to dispatch to serial
################################################################################
skel = SkeletonTriangulation(dpanel_model)

skel_panel_ids = get_panel_ids(skel)
_skel_panel_ids = get_skel_panel_ids(skel_panel_ids)

skel_geo_map = geo_map_func(_skel_panel_ids)
writevtk(skel,dir*"/distributed_skel_trian",append=false,geo_map=skel_geo_map)
