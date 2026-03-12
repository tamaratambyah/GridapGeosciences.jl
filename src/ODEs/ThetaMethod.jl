
##################
# Nonlinear case #
##################
function Gridap.ODEs.allocate_odecache(
  odeslvr::ThetaMethod, odeop::DAEODEOperator,
  t0::Real, us0::NTuple{1,AbstractVector}
)
  u0 = us0[1]
  us0N = (u0, u0)
  odeopcache = Gridap.ODEs.allocate_odeopcache(odeop, t0, us0N)

  uθ = copy(u0)

  sysslvrcache = nothing
  odeslvrcache = (uθ, sysslvrcache)

  dae_cache = nothing

  (odeslvrcache, odeopcache,dae_cache)
end

function Gridap.ODEs.ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::ThetaMethod, odeop::ODEOperator,
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache
)
println("DAE THETA")
  # Unpack inputs
  u0 = state0[1]
  odeslvrcache, odeopcache, dae_cache = odecache
  uθ, sysslvrcache = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr
  dt, θ = odeslvr.dt, odeslvr.θ

  dae_nl = odeop.dae_nl
  diagnostics = get_free_dof_values(odeopcache.diagnostics)
  daeop = DAENonlinearOperator(odeop,odeopcache,zero(u0),t0)

  # solve for initial diagnostics, and store in the ODE operator cache
  update!(daeop,u0,t0)
  dae_cache = solve!(diagnostics,dae_nl,daeop,dae_cache)
  Gridap.ODEs.update_odeopcache!(odeopcache, odeop, t0, diagnostics,diagnostics)

  # Define scheme
  x = stateF[1]
  dtθ = θ * dt
  tx = t0 + dtθ
  function usx(x)
    copy!(uθ, u0)
    axpy!(dtθ, x, uθ)
    (uθ, x)
  end
  ws = (dtθ, 1)

  # Update ODE operator cache
  Gridap.ODEs.update_odeopcache!(odeopcache, odeop, tx, diagnostics,diagnostics)

  # Create and solve stage operator
  stageop = Gridap.ODEs.NonlinearStageOperator(
    odeop, odeopcache,
    tx, usx, ws
  )

  sysslvrcache = Gridap.Algebra.solve!(x, sysslvr, stageop, sysslvrcache)

  # Update state
  tF = t0 + dt
  stateF = Gridap.ODEs._udate_theta!(stateF, state0, dt, x)

  # final diagnostic solve
  update!(daeop,stateF[1],tF)
  dae_cache = solve!(diagnostics,dae_nl,daeop,dae_cache)
  Gridap.ODEs.update_odeopcache!(odeopcache, odeop, t0 + dt, diagnostics)

  # Pack outputs
  odeslvrcache = (uθ, sysslvrcache)
  odecache = (odeslvrcache, odeopcache, dae_cache)
  (tF, stateF, odecache)
end
