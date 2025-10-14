using DrWatson
using Gridap
using GridapGeosciences
using Plots, LaTeXStrings


include("WaveEquation.jl")
# include("shallow_water.jl")
include("NLShallowWater.jl")

ps = [1]
ζs=[0.0]
n_ref_lvls = 4
ls=LUSolver()

################################################################################
#### Serial convergence test
################################################################################
ranks = [true]
nprocs = 1
WaveEquation.wave_convergence_test(ranks,nprocs,ζs,n_ref_lvls,ps,ls,true)
NLShallowWater.nonlinear_shallow_water_convergence_test(ranks,nprocs,ζs,n_ref_lvls,ps,ls,true,true)

# williamson2_convergence_test(linear_shallow_water_errors,n_ref_lvls)
# williamson2_convergence_test(nonlinear_shallow_water_errors,n_ref_lvls)
#

# ### linearised shallow water with arbitary functions
# depth(XYZ) = 1.0 + 0.1*exp(-( XYZ[2]^2 + XYZ[3]^2 ) )
# velocity(XYZ) = VectorValue(-XYZ[2],XYZ[1],0.0)
# coriolis(XYZ) = 2.0

# h = panel_to_cartesian(depth)
# vecX = velocity
# vX = panel_to_cartesian(tangent_vec(vecX))
# f = panel_to_cartesian(coriolis)

# linear_shallow_water_convergence_test(n_ref_lvls,h,vX,f,true)


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
