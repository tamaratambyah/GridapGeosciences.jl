
using MPI 
using PartitionedArrays

include("../L2ProjectionTests.jl")


MPI.Init()
nprocs = MPI.Comm_size(MPI.COMM_WORLD)
ranks = distribute_with_mpi(LinearIndices((prod(nprocs),)))

n_ref_lvls = 4

## Distributed model: 2D
models = get_distributed_refined_models(ranks,nprocs,n_ref_lvls)
L2_projection(models; _i_am_main=i_am_main(ranks))

### P4test model: 2D
models = get_octree_refined_models(ranks,n_ref_lvls)
L2_projection(models; _i_am_main=i_am_main(ranks))

### P4test model: 3D
models = get_3D_octree_refined_models(ranks,n_ref_lvls-1)
L2_projection(models; _i_am_main=i_am_main(ranks))
