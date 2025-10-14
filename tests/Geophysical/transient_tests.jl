using DrWatson
using MPI, PartitionedArrays


# include("TransientWaveEquation.jl")
include("TransientShallowWater.jl")

# options_cg = """
# -ksp_type cg
# -ksp_rtol 1.0e-6
# -ksp_converged_reason
# -ksp_monitor
# """
# with_debug() do distribute
#   TransientWaveEquation.main_transient(distribute,1;n_ref_lvls=3,options=options_cg,return_vtk=true)
# end

with_debug() do distribute
  TransientShallowWater.main_transient(distribute,6;n_ref_lvls=3,return_vtk=true)
end
