using MPI, PartitionedArrays
using GridapGeosciences
include("../AmbientSurfaceArea.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
ranks = distribute_with_mpi(LinearIndices((prod(nprocs),)))

n_ref_lvls = 4
radius = 1.0
models = get_distributed_ambient_refined_models(ranks,nprocs,n_ref_lvls,radius)
AmbientSurfaceArea.main(models;_i_am_main=i_am_main(ranks))
