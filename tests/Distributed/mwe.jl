using DrWatson
using Gridap
using GridapDistributed
using GridapP4est
using GridapGeosciences
using GridapGeosciences.Distributed

using MPI
using PartitionedArrays



MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

num_horizontal_uniform_refinements = 2
num_vertical_uniform_refinements = 1
omodel = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
	                                       num_horizontal_uniform_refinements=num_horizontal_uniform_refinements,
                                           num_vertical_uniform_refinements=num_vertical_uniform_refinements);

# omodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=num_horizontal_uniform_refinements)

panel_model = omodel.parametric_dmodel
panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(panel_model)

map(panel_ids,local_views(panel_model),local_views(Ω_panel)) do pid, lmodel, trian
  local_pids = get_panel_ids(lmodel)
  i_am_main(ranks) && println("Model panel ids: ", pid)
  i_am_main(ranks) && println("Local panel ids: ", pid)
  i_am_main(ranks) && println("length model pids = ", length(pid), "; num cells trian = ", num_cells(trian))
  i_am_main(ranks) && println("length local pids = ", length(local_pids), "; num cells local trian = ", num_cells(Triangulation(lmodel)))
end

### Broken for > 1 proc
metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
