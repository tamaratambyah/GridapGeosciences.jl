using DrWatson
using Gridap
using GridapGeosciences
using Plots, LaTeXStrings

include("../convergence_tools.jl")
include("TransientAdvectionSUPG.jl")
# include("AdvectionDGUpwinding.jl")
include("../Advection/advection_funcs.jl")

vX = vt #panel_to_cartesian(tangent_vec(vecX))
u = panel_to_cartesian(u0)
uvX = panel_to_cartesian(u0vecX)

n_ref_lvls = 4
CFL = 0.1
ps = [1]#[1,2,3]
tF=2*π
################################################################################
#### Serial convergence test
################################################################################
ls = LUSolver()
ranks = [true]
nprocs = 1

## SUPG
TransientAdvectionSUPG.transient_advection_supg_convergence_test(ranks,nprocs,u,vX,n_ref_lvls,ps,CFL,ls,tF,true)

## DG
# AdvectionDGUpwinding.advection_dg_convergence_test(ranks,nprocs,u,vX,uvX,n_ref_lvls,ps,ls,true)


################################################################################
#### Distributed convergence test
################################################################################
using GridapDistributed
using PartitionedArrays
using MPIPreferences
MPIPreferences.use_jll_binary()
using GridapSolvers, GridapPETSc


nprocs = 6

ranks = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end

ls = LUSolver()
## SUPG
TransientAdvectionSUPG.transient_advection_supg_convergence_test(ranks,nprocs,u,vX,n_ref_lvls,ps,CFL,ls,tF,true)

# ## DG
# ls = LUSolver()
# # using GridapSolvers
# # ls = GMRESSolver(10;Pr=JacobiLinearSolver(),maxiter=2000,verbose=i_am_main(ranks))
# AdvectionDGUpwinding.advection_dg_convergence_test(ranks,nprocs,u,vX,uvX,n_ref_lvls,ps,ls,true)
