"""
To run Galewsky test case
"""

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

include("../Geophysical/Galewsky.jl")
include("TransientThermalShallowWater.jl")



################################################################################
#### Main run for transient solution
################################################################################
function main_transient(distribute,nprocs;
  restart=false,n_ref_lvls=5,p_fe=1,CFL=0.1,return_vtk=true,_ε=1e-4,_soft=true)

  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("transient_tsw_Gakewsky_equation")

  vX = panel_to_cartesian(tangent_vec(u₀))
  f = panel_to_cartesian(f₀)
  h = panel_to_cartesian(h₀)
  B = panel_to_cartesian(B₀)

  ls_diag = CGSolver(JacobiLinearSolver();rtol=1-14,atol=1e-16,verbose=i_am_main(ranks),name="diagnostic_solver")
  ls_ode = CGSolver(JacobiLinearSolver();rtol=1-14,atol=1e-16,verbose=i_am_main(ranks),name="ode_solver")
  lss = (ls_ode,ls_diag)

  omodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=n_ref_lvls)
  panel_model = omodel.parametric_dmodel

  _dir = datadir("TransientThermalShallowWater_Galewsky")
  (i_am_main(ranks) && !isdir(_dir)) && mkdir(_dir)

  dir = _dir*"/sol_p$(p_fe)_nref$n_ref_lvls"
  (i_am_main(ranks) && !isdir(dir) && return_vtk) && mkdir(dir)

  transient_tsw_solver(panel_model,p_fe,dir,h,vX,f,B,_ε,_soft,CFL,lss,restart)
  post_process(panel_model,p_fe,dir,f,return_vtk)

  i_am_main(ranks) && println("--DONE--")
  @test true
end



function main_visualise(distribute,nprocs;
  n_ref_lvls=5,p_fe=1,return_vtk=true)

  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("Galewsky_visualise")

  f = panel_to_cartesian(f₀)

  omodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=n_ref_lvls)
  panel_model = omodel.parametric_dmodel

  _dir = datadir("TransientThermalShallowWater_Galewsky")

  dir = _dir*"/sol_p$(p_fe)_nref$n_ref_lvls"

  post_process(panel_model,p_fe,dir,f,return_vtk)

  i_am_main(ranks) && println("--DONE--")
  @test true
end
