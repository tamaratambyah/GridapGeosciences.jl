

##################
# Nonlinear case #
##################
function Gridap.ODEs.allocate_odecache(
  odeslvr::DIMRungeKutta, odeop::DAEODEOperator{<:AbstractQuasilinearODE},
  t0::Real, us0::NTuple{1,AbstractVector}
)
  println("DAE DIM RK solver")

  u0 = us0[1]
  us0N = (u0, u0)
  odeopcache = Gridap.ODEs.allocate_odeopcache(odeop, t0, us0N)

  ui_pre, ui = zero(u0), zero(u0)
  num_stages = length(get_nodes(odeslvr.tableau))
  slopes = [zero(u0) for _ in 1:num_stages]

  odeoptype = ODEOperatorType(odeop)
  has_explicit = odeoptype <: AbstractQuasilinearODE
  is_semilinear = odeoptype <: AbstractSemilinearODE
  mass_constant = Gridap.ODEs.is_form_constant(odeop, 1)
  reuse = (is_semilinear && mass_constant)

  J, r = nothing, nothing
  if has_explicit
    # Allocate J, r if there are explicit stages
    A = get_matrix(odeslvr.tableau)
    if any(i -> iszero(A[i, i]), axes(A, 2))
      J = Gridap.ODEs.allocate_jacobian(odeop, t0, us0N, odeopcache)
      r = Gridap.ODEs.allocate_residual(odeop, t0, us0N, odeopcache)
    end
  end

  sysslvrcaches = (nothing, nothing)
  odeslvrcache = (reuse, has_explicit, ui_pre, ui, slopes, J, r, sysslvrcaches)

  dae_cache = nothing

  (odeslvrcache, odeopcache, dae_cache)
end

function Gridap.ODEs.ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::DIMRungeKutta, odeop::DAEODEOperator{<:AbstractQuasilinearODE},
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache
)
  println("DAE DIM RK solver")

  # Unpack inputs
  u0 = state0[1]
  odeslvrcache, odeopcache, dae_cache = odecache
  reuse, has_explicit, ui_pre, ui, slopes, J, r, sysslvrcaches = odeslvrcache
  sysslvrcache_nl, sysslvrcache_l = sysslvrcaches

  # Unpack solver
  sysslvr_nl, sysslvr_l = odeslvr.sysslvr_nl, odeslvr.sysslvr_l
  dt, tableau = odeslvr.dt, odeslvr.tableau
  A, b, c = get_matrix(tableau), get_weights(tableau), get_nodes(tableau)

  dae_nl = odeop.dae_nl
  diagnostics = get_free_dof_values(odeopcache.diagnostics)
  daeop = DAENonlinearOperator(odeop,odeopcache,zero(u0),t0)

  # solve for initial diagnostics, and store in the ODE operator cache
  update!(daeop,u0,t0)
  dae_cache = solve!(diagnostics,dae_nl,daeop,dae_cache)
  Gridap.ODEs.update_odeopcache!(odeopcache, odeop, t0, diagnostics,diagnostics)

  for i in eachindex(c)
    # Define scheme
    x = slopes[i]
    tx = t0 + c[i] * dt
    copy!(ui_pre, u0)
    for j in 1:i-1
      axpy!(A[i, j] * dt, slopes[j], ui_pre)
    end

    # solve for stage diagnostics, and store in the ODE operator cache
    update!(daeop,ui_pre,tx)
    dae_cache = solve!(diagnostics,dae_nl,daeop,dae_cache)
    Gridap.ODEs.update_odeopcache!(odeopcache, odeop, tx, diagnostics)

    # Update ODE operator cache
    update_odeopcache!(odeopcache, odeop, tx)

    aii = A[i, i]

    # Define scheme
    function usx(x)
      copy!(ui, ui_pre)
      axpy!(aii * dt, x, ui)
      (ui, x)
    end
    ws = (aii * dt, 1)

    # Create and solve stage operator
    stageop = NonlinearStageOperator(
      odeop, odeopcache,
      tx, usx, ws
    )

    sysslvrcache_nl = Gridap.Algebra.solve!(x, sysslvr_nl, stageop, sysslvrcache_nl)
  end

  # Update state
  tF = t0 + dt
  stateF = Gridap.ODEs._update_dimrk!(stateF, state0, dt, slopes, b)

  # final diagnostic solve
  update!(daeop,stateF[1],tF)
  dae_cache = solve!(diagnostics,dae_nl,daeop,dae_cache)
  Gridap.ODEs.update_odeopcache!(odeopcache, odeop, t0 + dt, diagnostics)


  # Pack outputs
  sysslvrcaches = (sysslvrcache_nl, sysslvrcache_l)
  odeslvrcache = (reuse, has_explicit, ui_pre, ui, slopes, J, r, sysslvrcaches)
  odecache = (odeslvrcache, odeopcache, dae_cache)
  (tF, stateF, odecache)
end
