using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using MPI
using PartitionedArrays
using MPIPreferences
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra
using GridapGeosciences
using GridapPETSc
using Test

include("ThermogeostrophicBalanceTest.jl")
include("TransientThermalShallowWater.jl")

################################################################################
#### Main run for transient solution
################################################################################
function main_transient(distribute,nprocs;
  restart=false,options="",n_ref_lvls=4,p_fe=1,CFL=0.1,ζ=0.0,return_vtk=true,_ε=1e-4,_soft=false)

  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("transient_tsw_equation")

  h = panel_to_cartesian(h₀(ζ))
  vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
  f = panel_to_cartesian(f₀(ζ))
  B = panel_to_cartesian(B₀(ζ))

  ls_diag = CGSolver(JacobiLinearSolver();maxiter=1000,rtol=1-16,atol=1e-16,verbose=false,name="diagnostic_solver")
  ls_ode = GMRESSolver(10;Pr=JacobiLinearSolver(),maxiter=1000,rtol=1-14,verbose=i_am_main(ranks),name="ode_solver")
  lss = (ls_ode,ls_diag)

  radius = 1.0
  omodel = CubedSphere2DParametricOctreeDistributedDiscreteModel(ranks, radius; num_initial_uniform_refinements=n_ref_lvls)
  panel_model = omodel.parametric_dmodel

  _dir = datadir("TransientThermalShallowWater_checkpointing")
  (i_am_main(ranks) && !isdir(_dir)) && mkdir(_dir)

  dir = _dir*"/sol_p$(p_fe)_nref$n_ref_lvls"
  (i_am_main(ranks) && !isdir(dir) && return_vtk) && mkdir(dir)

  transient_tsw_solver(panel_model,p_fe,dir,h,vX,f,B,_ε,_soft,CFL,lss,restart)
  convergence_post_process(panel_model,p_fe,dir)
  post_process(panel_model,p_fe,dir,f,return_vtk)

  i_am_main(ranks) && println("--DONE--")
  @test true
end
