using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed

## should return ghost+owned panel ids
MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

omodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=1)
panel_model = omodel.parametric_dmodel
Ω_panel = Triangulation(panel_model)
get_panel_ids(panel_model)
get_panel_ids(Ω_panel)

cell_geo_map = geo_map_func(Ω_panel)
writevtk(Ω_panel,"test",append=false, geo_map=cell_geo_map)


ref_coarse_flags=map(ranks,partition(get_cell_gids(omodel.octree_dmodel))) do rank,indices
	flags=zeros(Cint,length(indices))
	flags.=refine_flag
	flags
end

omodel_adapted = Gridap.Adaptivity.adapt(omodel,ref_coarse_flags)
panel_model_adapted = omodel_adapted.parametric_dmodel
Ω_panel_adapted = Triangulation(panel_model_adapted)

get_panel_ids(panel_model_adapted)
get_panel_ids(Ω_panel_adapted)

cell_geo_map = geo_map_func(get_panel_ids(Ω_panel_adapted))
writevtk(Ω_panel_adapted,"test_adapted",append=false, geo_map=cell_geo_map)

using Gridap.Adaptivity
model = panel_model_adapted.models.item
get_panel_ids(model) ### incorrect

glue = get_adaptivity_glue(model)
pids = get_panel_ids(model.parent)
Gridap.Adaptivity.o2n_reindex(pids,glue)
