using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using Gridap.Algebra
using GridapDistributed

using DrWatson

include("../missing_overloads.jl")



MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

include("../convergence_tools.jl")
dir = datadir("Laplace3D")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)


n_ref_lvl = 4

p_fe = 1
ls = LUSolver()
return_vtk = true
fXYZ(XYZ) =  XYZ[1]*XYZ[2]*XYZ[3]
f = panel_to_cartesian(fXYZ)

for v_lvl in n_ref_v:-1:1
  models = get_3D_octree_refined_models(ranks,n_ref_lvl)
  errors,ns,dxs,slopes = h_convergence_test(models,laplace_beltrami_solver_3D,p_fe,dir,f,ls,return_vtk)
end
