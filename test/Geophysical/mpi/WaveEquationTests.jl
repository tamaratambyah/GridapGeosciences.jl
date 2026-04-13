using MPI, PartitionedArrays
using GridapGeosciences
include("../WaveEquation.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
ranks = distribute_with_mpi(LinearIndices((prod(nprocs),)))

n_ref_lvls = 4

## Distributed model: 2D
models = get_distributed_refined_models(ranks,nprocs,n_ref_lvls)
WaveEquationTests.main(models;_i_am_main=i_am_main(ranks))

### P4test model: 2D
models = get_octree_refined_models(ranks,n_ref_lvls)
WaveEquationTests.main(models;_i_am_main=i_am_main(ranks))

### P4test model: 3D
models = get_3D_octree_refined_models(ranks,n_ref_lvls-1)
WaveEquationTests.main(models;ps=[1],_i_am_main=i_am_main(ranks))


# with_mpi() do distribute
#   WaveEquationTests.main(distribute,nprocs)
# end
