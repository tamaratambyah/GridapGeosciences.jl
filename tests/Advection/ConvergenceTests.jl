using DrWatson
using Gridap
using GridapGeosciences
using Plots, LaTeXStrings

include("../convergence_tools.jl")
include("AdvectionSUPG.jl")
include("AdvectionDGUpwinding.jl")
include("advection_funcs.jl")


vX = panel_to_cartesian(tangent_vec(vecX))
u = panel_to_cartesian(u0)
uvX = panel_to_cartesian(u0vecX)

n_ref_lvls = 4
CFL = 0.1
ps = [1]#[1,2,3]

################################################################################
#### Serial convergence test
################################################################################
dir = datadir("SerialAdvectionTests")
!isdir(dir) && mkdir(dir)
ls = LUSolver()
ranks = [true]
nprocs = 1

## SUPG
AdvectionSUPG.advection_supg_convergence_test(ranks,nprocs,dir,u,vX,n_ref_lvls,ps,CFL,ls)

## DG
AdvectionDGUpwinding.advection_dg_convergence_test(ranks,nprocs,dir,u,vX,uvX,n_ref_lvls,ps,ls)




# transient_advection_supg_convergence_test(n_ref_lvls,u,vX,CFL,false) ### BROKEN

# transient_advection_dg_convergence_test(n_ref_lvls,u,vX,CFL,false) ### BROKEN



################################################################################
#### Distributed convergence test
################################################################################
using GridapDistributed
using PartitionedArrays
using MPIPreferences
MPIPreferences.use_jll_binary()
using GridapSolvers, GridapPETSc


nprocs = 6
dir = datadir("DistributedAdvectionTests")

ranks = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end
(i_am_main(ranks) && !isdir(dir) ) && mkdir(dir)

## SUPG
AdvectionSUPG.advection_supg_convergence_test(ranks,nprocs,dir,u,vX,n_ref_lvls,ps,CFL,ls)

## DG
AdvectionDGUpwinding.advection_dg_convergence_test(ranks,nprocs,dir,u,vX,uvX,n_ref_lvls,ps,ls)
