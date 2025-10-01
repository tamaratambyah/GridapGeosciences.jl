using DrWatson
using Gridap
using GridapDistributed
using PartitionedArrays
using MPIPreferences

using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers
MPIPreferences.use_jll_binary()

using GridapGeosciences
using GridapPETSc
include("convergence_tools.jl")

nprocs = 48
ranks = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end

### return distributed version of fine serial model
function get_distributed_panel_model(ranks,nprocs,n_ref_lvls::Int)

  # get refined models in serial
  s_models  = get_refined_models(n_ref_lvls)

  # s_panel_model = coarse_parametric_model()
  # for n in n_ref_lvls:-1:1
  #   s_panel_model = Gridap.Adaptivity.refine(s_panel_model)
  # end


  s_model_fine = s_models[1]
  spanel_ids = get_panel_ids(s_model_fine)

  part_to_cells = [PartitionedArrays.local_range(rank,nprocs,num_cells(s_model_fine)) for rank in 1:nprocs]

  # get the partition
  fine_cell_to_part = zeros(Int32,num_cells(s_model_fine))
  for (rank, cells) in enumerate(part_to_cells)
    fine_cell_to_part[cells] .= rank
  end

  # distribute the model
  dmodel = DiscreteModel(ranks,Adaptivity.get_model(s_model_fine),fine_cell_to_part)
  dpanel_ids, = distributed_panel_ids(dmodel,spanel_ids)
  panel_model = DistributedParametricDiscreteModel(dmodel,dpanel_ids)

 return panel_model
end


panel_model = get_distributed_panel_model(ranks,nprocs,6)

cell_geo_map = geo_map_func(get_panel_ids(panel_model))

dir = datadir("2D_CubedSphereRefactor")

writevtk(Triangulation(panel_model),dir*"/ambient_model_ref", append=false, compress=false,geo_map=cell_geo_map)

nprocs_per_panel = nprocs/6
nc_per_panel = sqrt(nc(panel_model))
nc_per_proc = nc_per_panel/nprocs_per_panel

nc_per_panel*nc_per_panel*6 == num_cells(panel_model)
