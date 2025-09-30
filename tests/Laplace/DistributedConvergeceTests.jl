using DrWatson
using Gridap
using GridapDistributed
using PartitionedArrays
using MPIPreferences

using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers
MPIPreferences.use_jll_binary()

using GridapGeosciences
using GridapSolvers
using Plots, LaTeXStrings
using GridapPETSc


include("../Distributed/convergence_tools.jl")
include("Helmholtz.jl")
include("LaplaceBeltrami.jl")
include("analytic_funcs.jl")



function main(convergence_func;nprocs,
  n_ref_lvls=4,
  analytic_funcs=analytic_funcs,
  ps=[1,2,3],
  out_loc="DistributedLaplaceTests",
  options=options_gmres,
  return_vtk=false)

  ranks = with_debug() do distribute
    distribute(LinearIndices((nprocs,)))
  end

  dir = datadir(out_loc)

  (i_am_main(ranks) && (return_vtk && !isdir(dir)) ) && mkdir(dir)

  GridapPETSc.Init(args=split(options))

  # ls = PETScLinearSolver(petsc_ls_from_options_g)
  ls = PETScLinearSolver(petsc_gmres_jacobi)
  # ls = GMRESSolver(10,Pr=JacobiLinearSolver(),rtol=1.e-12,verbose=true)
  # ls = LUSolver()

  convergence_func(ranks,nprocs,dir,analytic_funcs,n_ref_lvls,ps,ls,return_vtk)

  GridapPETSc.Finalize()
  GridapPETSc.gridap_petsc_gc()

  i_am_main(ranks) && println("--DONE--")

end


## Distributed Laplace Beltrami
for dict in [analytic_funcs, mapped_funcs, williamson_funcs]

  main(laplace_beltrami_convergence_test;nprocs=6,n_ref_lvls=4,
  analytic_funcs=dict,ps=[1,2,3],
  out_loc="DistributedLaplaceTests",options=options_gmres,return_vtk=false)

end


## Distributed Helmholtz
for dict in [analytic_funcs, mapped_funcs, williamson_funcs]

  main(helmholtz_convergence_test;nprocs=6,n_ref_lvls=4,
  analytic_funcs=dict,ps=[1,2,3],
  out_loc="DistributedLaplaceTests",options=options_gmres,return_vtk=false)

end
