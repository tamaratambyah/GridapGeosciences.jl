using MPI, PartitionedArrays
using GridapGeosciences
include("../Helmholtz.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
ranks = distribute_with_mpi(LinearIndices((prod(nprocs),)))

n_ref_lvls = 4

## Distributed model: 2D
models = get_distributed_refined_models(ranks,nprocs,n_ref_lvls)
HelmholtzTests.main(models;_i_am_main=i_am_main(ranks))

### P4test model: 2D
models = get_octree_refined_models(ranks,n_ref_lvls)
HelmholtzTests.main(models;_i_am_main=i_am_main(ranks))
