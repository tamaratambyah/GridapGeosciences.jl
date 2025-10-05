using DrWatson
using Gridap
using GridapGeosciences
using Plots, LaTeXStrings

include("../convergence_tools.jl")
include("Helmholtz.jl")
include("LaplaceBeltrami.jl")
include("MixedHelmholtz.jl")
include("analytic_funcs.jl")

n_ref_lvls = 4
ps = [1,2,3]

################################################################################
#### Serial convergence test
################################################################################
dir = datadir("SerialLaplaceTests")
!isdir(dir) && mkdir(dir)
ranks = [true]
nprocs = 1
for dict in [analytic_funcs]#, mapped_funcs, williamson_funcs]
  ### Laplace-Beltrami - requires zeromean functions
  LaplaceBeltrami.laplace_beltrami_convergence_test(ranks,nprocs,dir,dict,n_ref_lvls,ps)

  ### Helmholtz: primial formulation
  Helmholtz.helmholtz_convergence_test(ranks,nprocs,dir,dict,n_ref_lvls,ps)

  ### Helmholtz: mixed problem
  MixedHelmholtz.mixed_helmholtz_convergence_test(ranks,nprocs,dir,dict,n_ref_lvls,[1,2])
end




################################################################################
#### Distributed convergence test
################################################################################
using GridapDistributed
using PartitionedArrays
using MPIPreferences
MPIPreferences.use_jll_binary()
using GridapSolvers, GridapPETSc


nprocs = 6
dir = datadir("DistributedLaplaceTests")
cg = CGSolver(JacobiLinearSolver();maxiter=2000,verbose=true)
minres = MINRESSolver(;Pl=JacobiLinearSolver(),maxiter=5000,verbose=true)
ranks = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end

(i_am_main(ranks) && !isdir(dir) ) && mkdir(dir)


for dict in [analytic_funcs]#, mapped_funcs, williamson_funcs]

  LaplaceBeltrami.laplace_beltrami_convergence_test(ranks,nprocs,dir,dict,n_ref_lvls,ps,cg)

  Helmholtz.helmholtz_convergence_test(ranks,nprocs,dir,dict,n_ref_lvls,ps,cg)

  MixedHelmholtz.mixed_helmholtz_convergence_test(ranks,nprocs,dir,dict,n_ref_lvls,[1,2],minres)

end
