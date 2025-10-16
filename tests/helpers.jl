l2(e,dΩ) = sum(∫( e⋅e )dΩ)
nc(panel_model) = num_cells(panel_model)/6 ## nc = num cells per panel
dx(nc) = sqrt( 4*π*RADIUS^2 / (6*sqrt(nc)^2) )
nref(nc) = Int(log2(sqrt(nc))) ## level of refinement

using Printf
using GridapSolvers
using GridapSolvers.SolverInterfaces: ConvergenceLog

function GridapSolvers.SolverInterfaces.init!(log::ConvergenceLog{T},r0::T) where T
  GridapSolvers.SolverInterfaces.reset!(log)
  log.residuals[1] = r0
  # if log.verbose > SOLVER_VERBOSE_LOW
    header =  " Starting $(log.name) solver "
    println(GridapSolvers.SolverInterfaces.get_tabulation(log,0),rpad(string(repeat('-',15),header),55,'-'))
    t = GridapSolvers.SolverInterfaces.get_tabulation(log)
    msg = @sprintf("> Iteration %3i - Residuals: %.2e,   %.2e ", 0, r0, 1)
    println(t,msg)
  # end
  return GridapSolvers.SolverInterfaces.finished(log.tols,log.num_iters,r0,1.0)
end

function GridapSolvers.SolverInterfaces.finalize!(log::ConvergenceLog{T},r::T) where T
  r_rel = r / log.residuals[1]
  flag  = GridapSolvers.SolverInterfaces.finished_flag(log.tols,log.num_iters,r,r_rel)
  # if log.verbose > SOLVER_VERBOSE_NONE
    t = GridapSolvers.SolverInterfaces.get_tabulation(log,0)
    println(t,"Solver $(log.name) finished with reason $(flag)")
    msg = @sprintf("Iterations: %3i - Residuals: %.2e,   %.2e ", log.num_iters, r, r_rel)
    println(t,msg)
    # if log.verbose > SOLVER_VERBOSE_LOW
      footer = " Exiting $(log.name) solver "
      println(t,rpad(string(repeat('-',15),footer),55,'-'))
    # end
  # end
  return flag
end
