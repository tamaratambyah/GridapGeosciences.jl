import Gridap.ODEs: EXRungeKutta
import Gridap.ODEs: ODEOperator
import Gridap.ODEs: AbstractTableau
import Gridap.ODEs: AbstractQuasilinearODE
using LinearAlgebra
import Gridap.ODEs: get_weights, get_nodes,get_matrix
import Gridap.ODEs: RungeKutta
import Gridap.Algebra: NonlinearSolver
import Gridap.ODEs: ODESolver
import Gridap.ODEs: TableauType, ExplicitTableau
import Gridap.ODEs: SemilinearODE

###############
# Linear case #
###############
function Gridap.ODEs.allocate_odecache(
  odeslvr::EXRungeKutta, odeop::DAEODEOperator{<:AbstractQuasilinearODE},
  t0::Real, us0::NTuple{1,AbstractVector}
)
  u0 = us0[1]
  us0N = (u0, u0)
  odeopcache = Gridap.ODEs.allocate_odeopcache(odeop, t0, us0N)

  ui_pre = zero(u0)
  num_stages = length(get_nodes(odeslvr.tableau))
  slopes = [zero(u0) for _ in 1:num_stages]

  is_semilinear = ODEOperatorType(odeop) <: AbstractSemilinearODE
  constant_mass = Gridap.ODEs.is_form_constant(odeop, 1)
  reuse = (is_semilinear && constant_mass)

  J = Gridap.ODEs.allocate_jacobian(odeop, t0, us0N, odeopcache)
  r = Gridap.ODEs.allocate_residual(odeop, t0, us0N, odeopcache)

  sysslvrcache = nothing
  odeslvrcache = (reuse, ui_pre, slopes, J, r, sysslvrcache)

  dae_cache = nothing

  (odeslvrcache, odeopcache, dae_cache)
end

# function LinearAlgebra.axpy!(α,x::PVector,y::PVector)
#   map(PartitionedArrays.partition(x),PartitionedArrays.partition(y)) do x,y
#     LinearAlgebra.axpy!(α,x,y)
#   end
#   consistent!(y) |> fetch
#   return y
# end

# function Base.copy!(x::PVector,y::PVector)
#   map(PartitionedArrays.own_values(x),PartitionedArrays.own_values(y)) do x,y
#     copy!(x,y)
#   end
#   consistent!(x) |> fetch
#   return x
# end

function Gridap.ODEs.ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::EXRungeKutta, odeop::DAEODEOperator{<:AbstractQuasilinearODE},
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache
)

  println("new RK solver")

  # Unpack inputs
  u0 = state0[1]
  odeslvrcache, odeopcache, dae_cache = odecache
  reuse, ui_pre, slopes, J, r, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr
  dt, tableau = odeslvr.dt, odeslvr.tableau
  A, b, c = get_matrix(tableau), get_weights(tableau), get_nodes(tableau)


  dae_nl = odeop.dae_nl
  diagnostics = get_free_dof_values(odeopcache.diagnostics)
  daeop = DAENonlinearOperator(odeop,odeopcache,zero(u0),t0)


  for i in eachindex(c)

    tx = t0 + c[i] * dt
    ###################
    # Diagnostics     #
    ###################
    # make ui
    copy!(ui_pre, u0)
    for j in 1:i-1
      axpy!(A[i, j] * dt, slopes[j], ui_pre)
    end

    # diagnostics
    update!(daeop,ui_pre,tx)
    dae_cache = solve!(diagnostics,dae_nl,daeop,dae_cache)

    # Update ODE operator cache
    Gridap.ODEs.update_odeopcache!(odeopcache, odeop, tx, diagnostics)

    ###################
    # Runge Kutta     #
    ###################
    # Define scheme
    # Set x to zero to split jacobian and residual
    x = slopes[i]
    fill!(x, zero(eltype(x)))

     # make ui again - hack to avoid copy errors
    copy!(ui_pre, u0)
    for j in 1:i-1
      axpy!(A[i, j] * dt, slopes[j], ui_pre)
    end
    usx = (ui_pre, x)
    ws = (0, 1)


    # Create and solve stage operator
    stageop = LinearStageOperator(
      odeop, odeopcache,
      tx, usx, ws,
      J, r, reuse, sysslvrcache
    )

    sysslvrcache = Gridap.Algebra.solve!(x, sysslvr, stageop, sysslvrcache)
  end

  # Update state
  tF = t0 + dt
  stateF = Gridap.ODEs._update_exrk!(stateF, state0, dt, slopes, b)

  # final diagnostic solve
  update!(daeop,stateF[1],tF)
  dae_cache = solve!(diagnostics,dae_nl,daeop,dae_cache)
  Gridap.ODEs.update_odeopcache!(odeopcache, odeop, t0 + dt, diagnostics)


  # Pack outputs
  odeslvrcache = (reuse, ui_pre, slopes, J, r, sysslvrcache)
  odecache = (odeslvrcache, odeopcache, dae_cache)
  (tF, stateF, odecache)
end
