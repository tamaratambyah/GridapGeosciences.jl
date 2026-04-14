using MPI, PartitionedArrays
using GridapGeosciences
include("../HodgeLaplacian_vector.jl")

MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
ranks = distribute_with_mpi(LinearIndices((prod(nprocs),)))

n_ref_lvls = 3

### P4test model: 3D
models = get_3D_octree_refined_models(ranks,n_ref_lvls)
HodgeLaplacianVectorTests.main(models;_i_am_main=i_am_main(ranks))

# with_mpi() do distribute
#   HodgeLaplacianVectorTests.main(distribute,nprocs)
# end
