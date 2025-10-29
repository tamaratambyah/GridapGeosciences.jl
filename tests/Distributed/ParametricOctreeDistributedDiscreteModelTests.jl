using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed

## should return ghost+owned panel ids
function GridapGeosciences.get_panel_ids(omodel::ParametricOctreeDistributedDiscreteModel)
  panel_ids = map(local_views(omodel.parametric_dmodel)) do model
    return get_panel_ids(model)
  end
  return panel_ids
end

model   = coarse_parametric_model()
fmodel  = refine(model)

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

parametric_octree_dmodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=0)
get_panel_ids(parametric_octree_dmodel)
Triangulation(parametric_octree_dmodel)

map(local_views(parametric_octree_dmodel.parametric_dmodel)) do model
	panel_ids = get_panel_ids(model)
    cell_geo_map = lazy_map(p -> ForwardMap(p), panel_ids)
	writevtk(Triangulation(model),"test", geo_map=cell_geo_map)
end

ref_coarse_flags=map(ranks,partition(get_cell_gids(parametric_octree_dmodel.octree_dmodel))) do rank,indices
	flags=zeros(Cint,length(indices))
	flags.=refine_flag
	flags
end

parametric_octree_model_adapted = Gridap.Adaptivity.adapt(parametric_octree_dmodel,ref_coarse_flags)
get_panel_ids(parametric_octree_model_adapted)
get_panel_ids(parametric_octree_model_adapted.parametric_dmodel)
map(local_views(parametric_octree_model_adapted.parametric_dmodel)) do model
	panel_ids = get_panel_ids(model)
	println(typeof(model))
	cell_geo_map = lazy_map(p -> ForwardMap(p), panel_ids)
	writevtk(Triangulation(model),"test_adapted", geo_map=cell_geo_map)
end
