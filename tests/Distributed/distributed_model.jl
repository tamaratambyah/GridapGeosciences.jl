using DrWatson
using Gridap
using GridapDistributed
using PartitionedArrays
using MPIPreferences

using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers
MPIPreferences.use_jll_binary()

using GridapGeosciences

nprocs = 48
ranks = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end




panel_model = get_distributed_panel_model(ranks,nprocs,6)

cell_geo_map = geo_map_func(get_panel_ids(panel_model))

dir = datadir("2D_CubedSphereRefactor")

writevtk(Triangulation(panel_model),dir*"/ambient_model_ref", append=false, compress=false,geo_map=cell_geo_map)

nprocs_per_panel = nprocs/6
nc_per_panel = sqrt(num_cells(panel_model)/6)
nc_per_proc = nc_per_panel/nprocs_per_panel

nc_per_panel*nc_per_panel*6 == num_cells(panel_model)
