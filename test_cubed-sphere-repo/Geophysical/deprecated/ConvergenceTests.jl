using DrWatson
using Gridap
using GridapGeosciences
using Plots, LaTeXStrings


include("WaveEquation.jl")
include("LinearisedShallowWater.jl")
include("ShallowWater.jl")

ps = [1]
ζs = [0.0]
n_ref_lvls = 4
ls = LUSolver()
CFL = 0.1

################################################################################
#### Serial convergence test
################################################################################
ranks = [true]
nprocs = 1
WaveEquation.wave_convergence_test(ranks,nprocs,ζs,n_ref_lvls,ps,ls,CFL,true)
LinearisedShallowWater.linear_shallow_water_convergence_test(ranks,nprocs,ζs,n_ref_lvls,ps,ls,CFL,true,true)
ShallowWater.nonlinear_shallow_water_convergence_test(ranks,nprocs,ζs,n_ref_lvls,ps,ls,CFL,true,true)


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

WaveEquation.wave_convergence_test(ranks,nprocs,ζs,n_ref_lvls,ps,ls,true)

LinearisedShallowWater.linear_shallow_water_convergence_test(ranks,nprocs,ζs,n_ref_lvls,ps,ls,CFL,true,true)

ShallowWater.nonlinear_shallow_water_convergence_test(ranks,nprocs,ζs,n_ref_lvls,ps,ls,CFL,true,true)
