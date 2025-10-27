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

num_horizontal_uniform_refinements = 2
num_vertical_uniform_refinements = 0
a = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
	                                       num_horizontal_uniform_refinements=num_horizontal_uniform_refinements,
                                           num_vertical_uniform_refinements=num_vertical_uniform_refinements);

include("forward_map_3D.jl")

## Try plotting in distributed
model = a.parametric_dmodel
Ω_panel =  Triangulation(model)

cell_geo_map = geo_map_func_3D(Ω_panel)
writevtk(Ω_panel,dir*"/extruded_model",append=false,geo_map=cell_geo_map)
