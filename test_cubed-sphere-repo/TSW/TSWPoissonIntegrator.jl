"""
    struct PoissonIntegrator <: ODESolver end

"""
struct TSWPoissonIntegrator <: ODESolver
  sysslvr::NonlinearSolver
  diagslvr
  dt::Real
end


function Gridap.ODEs.allocate_odecache(
  odeslvr::TSWPoissonIntegrator, odeop::TSWODEOperator,
  t0::Real, us0::NTuple{1,AbstractVector},
  bhs, yhs, whs, zhs
)
  # println("allocate TSW ode cache")

  ts = (t0,odeslvr.dt)
  u0 = us0[1]
  us0N = (u0, u0)

  y = yhs[1]
  b = bhs[1]
  w = whs[1]
  z = zhs[1]

  odeopcache = Gridap.ODEs.allocate_odeopcache(odeop, t0, us0N, y,b,w,z)

  J_prog = allocate_prognostic_jacobian(odeop,ts,u0,odeopcache, u0, b, y, w, z, b)

  const_jac = is_const_jac(odeop)

  if const_jac
    println("The prognostic jacobian will be computed once per time step")
  else
    println("The prognostic will be recomputed at every newton iteration")
  end

  sysslvrcache = nothing
  odeslvrcache = (sysslvrcache, const_jac, J_prog)

  (odeslvrcache, odeopcache)
end


function Gridap.ODEs.ode_start(
  odeslvr::TSWPoissonIntegrator, odeop::TSWODEOperator,
  t0::Real, us0::Tuple{Vararg{AbstractVector}},
  odecache,
  bh, yh, wh, zh
)
  # println("tsw ode start")

  state0 = copy.(us0)
  diagnosticsb0 = copy.(bh)
  diagnostics0 = copy.(yh)
  diagonsticsw0 = copy.(wh)
  diagnosticsz0 = copy.(zh)

  # compute the prognostic jacobian
  u0 = state0[1]
  b0 =   diagnosticsb0[1]
  y0 = diagnostics0[1]
  w0 = diagonsticsw0[1]
  z0 = diagnosticsz0[1]

  odeslvrcache, odeopcache = odecache
  sysslvrcache, const_jac, J_prog = odeslvrcache
  ts = (t0,odeslvr.dt)
  prognostic_jacobian!(J_prog, odeop, ts, u0, odeopcache, u0, b0, y0, w0, z0, b0)

  # Pack outputs
  odeslvrcache = (sysslvrcache, const_jac, J_prog)
  odecache = (odeslvrcache, odeopcache)

  (state0, odecache, diagnosticsb0, diagnostics0,diagonsticsw0,diagnosticsz0)
end


function Gridap.ODEs.ode_march!(
  stateF::NTuple{1,AbstractVector},
  odeslvr::TSWPoissonIntegrator, odeop::TSWODEOperator,
  t0::Real, state0::NTuple{1,AbstractVector},
  odecache,
  diagnosticsbF,diagnosticsF,diagnosticswF,diagnosticszF,
  diagnosticsb0,diagnostics0,diagnosticsw0, diagnosticsz0
)
  # println("tsw ode march")


  # Unpack inputs
  u0 = state0[1]
  b0 = diagnosticsb0[1]
  y0 = diagnostics0[1]
  w0 = diagnosticsw0[1]
  z0 = diagnosticsz0[1]

  odeslvrcache, odeopcache = odecache
  sysslvrcache, const_jac, J_prog = odeslvrcache

  # Unpack solver
  sysslvr = odeslvr.sysslvr
  diagslvrs = odeslvr.diagslvr
  diagslvr1,diagslvr2 = diagslvrs

  dt = odeslvr.dt

  # Define scheme
  x = stateF[1]
  tx = t0
  ts = (tx,dt)

  y = diagnosticsF[1]
  b = diagnosticsbF[1]
  w = diagnosticswF[1]
  z = diagnosticszF[1]

  # Update ODE operator cache
  Gridap.ODEs.update_odeopcache!(odeopcache, odeop, tx)


  ############################################################################
  prognostic_jacobian!(J_prog, odeop, ts, u0, odeopcache, u0, b0, y0, w0, z0, b0)


  # # Create and solve stage operator
  stageop = TSWProgNonlinearOperator(odeop, odeopcache, ts, u0, diagslvrs, J_prog, const_jac,
                b0, y,w,z, b)

  sysslvrcache = Gridap.Algebra.solve!(x, sysslvr, stageop, sysslvrcache)



  # # Final solve for diagnostics
  diagslvrcache = odeopcache.diagslvrcache
  diagslvrcachew = odeopcache.diagslvrcachew
  diagslvrcachez = odeopcache.diagslvrcachez
  diagslvrcacheb = odeopcache.diagslvrcacheb


  diagbop = TSWDiagNonlinearOperatorb(odeop, odeopcache, ts, u0, x, b0)
  diagslvrcacheb = Gridap.Algebra.solve!(b, diagslvr2, diagbop, diagslvrcacheb)
  update_diagnosticsb!(odeopcache, b, diagslvrcacheb)


  diagop = TSWDiagNonlinearOperator(odeop, odeopcache, ts, u0, x)
  diagslvrcache = Gridap.Algebra.solve!(y, diagslvr1, diagop, diagslvrcache)
  update_diagnostics!(odeopcache, y, diagslvrcache)


  # # # solve for bhat
  # diagopw = TSWDiagNonlinearOperatorw(odeop, odeopcache, ts, b0, b, u0, x)
  # diagslvrcachew = Gridap.Algebra.solve!(w, diagslvr1, diagopw, diagslvrcachew)
  # update_diagnosticsw!(odeopcache, w, diagslvrcachew)

  # # solve for btilde
  diagopz = TSWDiagNonlinearOperatorz(odeop, odeopcache, ts, w, b0, b, y)
  diagslvrcachez = Gridap.Algebra.solve!(z, diagslvr1, diagopz, diagslvrcachez)
  update_diagnosticsz!(odeopcache, z, diagslvrcachez)

  copy_to_cache!(odeopcache, b,y,w,z)

  tF = t0 + dt
  # Update state
  stateF,diagnosticsbF,diagnosticsF,diagnosticswF,diagnosticszF = update_poisson!(stateF, x,
                          diagnosticsbF,diagnosticsF,diagnosticswF,diagnosticszF,
                            b,y,w,z)



  # Pack outputs
  odeslvrcache = (sysslvrcache, const_jac, J_prog)
  odecache = (odeslvrcache, odeopcache)
  (tF, stateF, odecache, diagnosticsbF,diagnosticsF,diagnosticswF,diagnosticszF)



end


function Gridap.ODEs.ode_finish!(
  uF::AbstractVector,
  odeslvr::TSWPoissonIntegrator, odeop::TSWODEOperator,
  t0::Real, tF, stateF::Tuple{Vararg{AbstractVector}},
  odecache,
  bF,yF,wF,zF,
  diagnosticsbF,diagnosticsF,diagnosticswF,diagnosticszF
  )
  # println("ODE finish")

  copy!(uF, first(stateF))

  copy!(bF, first(diagnosticsbF))
  copy!(yF, first(diagnosticsF) )
  copy!(wF, first(diagnosticswF) )
  copy!(zF, first(diagnosticszF))

  (uF, odecache, bF,yF,wF,zF)
end

#########
# Utils #
#########
function update_poisson!(
  stateF::NTuple{1,AbstractVector}, x::AbstractVector,
  diagnosticsbF,diagnosticsF,diagnosticswF,diagnosticszF,
  b,y,w,z
)
  uF = stateF[1]
  copy!(uF, x)

  bF = diagnosticsbF[1]
  copy!(bF,b)

  yF = diagnosticsF[1]
  copy!(yF,y)

  wF = diagnosticswF[1]
  copy!(wF,w)

  zF = diagnosticszF[1]
  copy!(zF,z)


  ((uF,),(bF,),(yF,),(wF,),(zF,))
end
