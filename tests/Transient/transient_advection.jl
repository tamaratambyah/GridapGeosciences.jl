using MPI, PartitionedArrays

include("TransientAdvectionSUPG.jl")
include("TransientWaveEquation.jl")
options_gmres = """
-ksp_type gmres
-ksp_rtol 1.0e-6
-ksp_converged_reason
-ksp_monitor
"""
with_debug() do distribute
  TransientAdvectionSUPG.main_transient(distribute;nprocs=6,options=options_gmres,
    n_ref_lvls=6,p_fe=1,CFL=0.1,tF=0.1,return_vtk=true)
end



options_cg = """
-ksp_type cg
-ksp_rtol 1.0e-6
-ksp_converged_reason
-ksp_monitor
"""
with_debug() do distribute
  TransientWaveEquation.main_transient(distribute,1;n_ref_lvls=3,options=options_cg,return_vtk=true)
end
