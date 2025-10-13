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
ps = [1]#[1,2,3]

################################################################################
#### Serial convergence test
################################################################################
ranks = [true]
nprocs = 1
ls = LUSolver()
for dict in [analytic_funcs]#, mapped_funcs, williamson_funcs]
  ### Laplace-Beltrami - requires zeromean functions
  LaplaceBeltrami.laplace_beltrami_convergence_test(ranks,nprocs,dict,n_ref_lvls,ps,ls,true)

  ### Helmholtz: primial formulation
  Helmholtz.helmholtz_convergence_test(ranks,nprocs,dict,n_ref_lvls,ps,ls,true)

  ### Helmholtz: mixed problem
  MixedHelmholtz.mixed_helmholtz_convergence_test(ranks,nprocs,dict,n_ref_lvls,[1,2],ls,true)
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
cg = CGSolver(JacobiLinearSolver();maxiter=2000,verbose=true)
minres = MINRESSolver(;Pl=JacobiLinearSolver(),maxiter=5000,verbose=true)
ranks = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end


for dict in [analytic_funcs]#, mapped_funcs, williamson_funcs]

  LaplaceBeltrami.laplace_beltrami_convergence_test(ranks,nprocs,dict,n_ref_lvls,ps,cg,true)

  Helmholtz.helmholtz_convergence_test(ranks,nprocs,dict,n_ref_lvls,ps,cg,true)

  MixedHelmholtz.mixed_helmholtz_convergence_test(ranks,nprocs,dict,n_ref_lvls,[1],minres,true)

end
