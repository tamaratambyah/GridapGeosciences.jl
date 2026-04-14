using MPI, PartitionedArrays
using GridapGeosciences
include("../Hierarchy.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
ranks = distribute_with_mpi(LinearIndices((prod(nprocs),)))

n_ref_lvls = 3

### P4test model: 2D
coarse_model = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=0)
HierarchyTest.main(coarse_model,n_ref_lvls)
