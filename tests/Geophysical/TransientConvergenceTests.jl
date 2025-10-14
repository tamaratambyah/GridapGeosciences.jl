using DrWatson
using Gridap
using GridapGeosciences
using Plots, LaTeXStrings


include("TransientShallowWater.jl")


ζs = [0.0]
n_ref_lvls = 3
ps = [1]
ls = LUSolver()
CFL = 0.1

################################################################################
#### Serial convergence test
################################################################################
ranks = [true]
nprocs = 1
TransientShallowWater.transient_shallow_water_convergence_test(ranks,nprocs,ζs,n_ref_lvls,ps,ls,CFL,true)



with_debug() do distribute
  main(distribute,1)
end
