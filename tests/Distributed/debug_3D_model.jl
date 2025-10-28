using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap

using DrWatson
dir = datadir("Distributed")
!isdir(dir) && mkdir(dir)

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

num_horizontal_uniform_refinements = 3
num_vertical_uniform_refinements = 1
o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
	                                       num_horizontal_uniform_refinements=num_horizontal_uniform_refinements,
                                           num_vertical_uniform_refinements=num_vertical_uniform_refinements);


panel_model = o3model.parametric_dmodel
Ω_panel =  Triangulation(panel_model)
panel_ids = get_panel_ids(Ω_panel)
panel_ids = get_panel_ids(panel_model)

## Try plotting in distributed
cell_geo_map = geo_map_func(Ω_panel)
writevtk(Ω_panel,dir*"/extruded_model",append=false,geo_map=cell_geo_map)
