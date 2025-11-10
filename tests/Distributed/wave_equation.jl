using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using PartitionedArrays
using MPI
using GridapP4est

using GridapGeosciences
using GridapGeosciences.Distributed

using Test


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

include("../missing_overloads.jl")
include("../convergence_tools.jl")
include("williamson_funcs_3D.jl")




dir = datadir("Wave3D")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

n_ref_lvl = 4

p_fe = 1
ls = LUSolver()
return_vtk = true
ζ = 0.0
h = panel_to_cartesian(h₀(ζ))
vX = panel_to_cartesian(tangent_vec(u₀(ζ)))

models = get_3D_octree_refined_models(ranks,n_ref_lvl)
errors,ns,dxs,slopes = h_convergence_test(models,wave_solver_3D,p_fe,dir,h,vX,ls,return_vtk)
