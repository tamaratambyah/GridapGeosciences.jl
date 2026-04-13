
using MPI, PartitionedArrays
using GridapGeosciences
include("../TransientShallowWater.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
ranks = distribute_with_mpi(LinearIndices((prod(nprocs),)))

## Distributed model: 2D
models = get_distributed_refined_models(ranks,nprocs,3)
TransientShallowWaterTests.main(models[1];_i_am_main=i_am_main(ranks))

### P4test model: 2D
omodel = ParametricOctreeDistributedDiscreteModel(ranks;
  num_initial_uniform_refinements=3)
panel_model = omodel.parametric_dmodel
TransientShallowWaterTests.main(panel_model;_i_am_main=i_am_main(ranks))
