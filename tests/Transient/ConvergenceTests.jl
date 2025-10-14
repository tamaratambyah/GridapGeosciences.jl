using DrWatson
using Gridap
using GridapGeosciences
using Plots, LaTeXStrings

include("../convergence_tools.jl")
include("TransientAdvectionSUPG.jl")
include("TransientAdvectionDGUpwinding.jl")
include("../Advection/advection_funcs.jl")

vX = panel_to_cartesian(tangent_vec(vecX))
u = panel_to_cartesian(u0)

n_ref_lvls = 4
CFL = 0.1
ps = [1]#[1,2,3]
tF = 2*π
################################################################################
#### Serial convergence test
################################################################################
ls = LUSolver()
ranks = [true]
nprocs = 1

## SUPG
TransientAdvectionSUPG.transient_advection_supg_convergence_test(ranks,nprocs,u,vt,n_ref_lvls,ps,CFL,ls,tF,true)

## DG
TransientAdvectionDGUpwinding.transient_advection_dg_convergence_test(ranks,nprocs,u,vX,n_ref_lvls,ps,CFL,ls,tF,true)


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
TransientAdvectionSUPG.transient_advection_supg_convergence_test(ranks,nprocs,u,vt,n_ref_lvls,ps,CFL,ls,tF,true)

# ## DG
TransientAdvectionDGUpwinding.transient_advection_dg_convergence_test(ranks,nprocs,u,vX,n_ref_lvls,ps,CFL,ls,tF,true)
