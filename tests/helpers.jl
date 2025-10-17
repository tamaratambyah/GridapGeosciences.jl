l2(e,dΩ) = sum(∫( e⋅e )dΩ)
nc(panel_model) = num_cells(panel_model)/6 ## nc = num cells per panel
dx(nc) = sqrt( 4*π*RADIUS^2 / (6*sqrt(nc)^2) )
nref(nc) = Int(log2(sqrt(nc))) ## level of refinement

using Printf
using GridapSolvers
using GridapSolvers.SolverInterfaces: ConvergenceLog


function GridapSolvers.SolverInterfaces.update!(log::ConvergenceLog{T},r::T) where T
  log.num_iters += 1
  log.residuals[log.num_iters+1] = r
  r_rel = r / log.residuals[1]
  if log.verbose > SOLVER_VERBOSE_LOW
    t = get_tabulation(log)
    # msg = @sprintf("> Iteration %3i - Residuals: %.2e,   %.2e ", log.num_iters, r, r_rel)
    # println(t,msg)
  end
  return finished(log.tols,log.num_iters,r,r_rel)
end
