
######################
# GenericODESolution #
######################
"""
    struct GenericODESolution <: ODESolution end

Generic wrapper for the evolution of an `ODEOperator` with an `ODESolver`.
"""
struct TSWGenericODESolution <: ODESolution
  odeslvr::ODESolver
  odeop::ODEOperator
  t0::Real
  tF::Real
  us0::Tuple{Vararg{AbstractVector}}
  bhs
  yhs
  whs
  zhs
end

function Base.iterate(odesltn::TSWGenericODESolution)
  odeslvr, odeop = odesltn.odeslvr, odesltn.odeop
  t0, us0 = odesltn.t0, odesltn.us0
  bhs, yhs, whs, zhs = odesltn.bhs, odesltn.yhs, odesltn.whs, odesltn.zhs

  # Allocate cache
  odecache = Gridap.ODEs.allocate_odecache(odeslvr, odeop, t0, us0, bhs, yhs, whs, zhs)

  # Starting procedure
  state0, odecache, diagnosticsb0, diagnostics0, diagonsticsw0, diagnosticsz0 = Gridap.ODEs.ode_start(
    odeslvr, odeop,
    t0, us0,
    odecache,
    bhs, yhs, whs, zhs
  )

  # Marching procedure
  stateF = copy.(state0)
  diagnosticsbF = copy.(diagnosticsb0)
  diagnosticsF = copy.(diagnostics0)
  diagnosticswF = copy.(diagonsticsw0)
  diagnosticszF = copy.(diagnosticsz0)

  tF, stateF, odecache, diagnosticsbF,diagnosticsF,diagnosticswF,diagnosticszF = Gridap.ODEs.ode_march!(
    stateF,
    odeslvr, odeop,
    t0, state0,
    odecache,
    diagnosticsbF,diagnosticsF,diagnosticswF,diagnosticszF,
    diagnosticsb0,diagnostics0,diagonsticsw0, diagnosticsz0
  )

  # Finishing procedure
  uF = copy(first(us0))
  bF = copy(first(bhs))
  yF = copy(first(yhs))
  wF = copy(first(whs))
  zF = copy(first(zhs))

  uF, odecache,bF,yF,wF,zF = Gridap.ODEs.ode_finish!(
    uF,
    odeslvr, odeop,
    t0, tF, stateF,
    odecache,
    bF,yF,wF,zF,
    diagnosticsbF,diagnosticsF,diagnosticswF,diagnosticszF
  )

  # Update iterator
  data = (tF, uF)
  state = (tF, stateF, state0, uF, odecache,
          bF,yF,wF,zF ,
          diagnosticsbF,diagnosticsF,diagnosticswF,diagnosticszF,
          diagnosticsb0, diagnostics0, diagonsticsw0, diagnosticsz0)
  (data, state)
end

function Base.iterate(odesltn::TSWGenericODESolution, state)
  odeslvr, odeop = odesltn.odeslvr, odesltn.odeop
  t0, state0, stateF, uF, odecache,
  bF,yF,wF,zF,
  diagnosticsb0, diagnostics0,diagonsticsw0,diagnosticsz0,
  diagnosticsbF,diagnosticsF,diagnosticswF,diagnosticszF  = state

  if t0 >= odesltn.tF - 100 * eps()
    return nothing
  end

  # Marching procedure
  tF, stateF, odecache, diagnosticsbF,diagnosticsF,diagnosticswF,diagnosticszF = Gridap.ODEs.ode_march!(
    stateF,
    odeslvr, odeop,
    t0, state0,
    odecache,
    diagnosticsbF,diagnosticsF,diagnosticswF,diagnosticszF,
    diagnosticsb0,diagnostics0,diagonsticsw0, diagnosticsz0
  )


  # Finishing procedure
  uF, odecache,bF,yF,wF,zF = Gridap.ODEs.ode_finish!(
    uF,
    odeslvr, odeop,
    t0, tF, stateF,
    odecache,
    bF,yF,wF,zF,
    diagnosticsbF,diagnosticsF,diagnosticswF,diagnosticszF
  )

  # Update iterator
  data = (tF, uF)
  state = (tF, stateF, state0, uF, odecache,
              bF,yF,wF,zF,
              diagnosticsbF,diagnosticsF,diagnosticswF,diagnosticszF,
              diagnosticsb0, diagnostics0,diagonsticsw0,diagnosticsz0)
  (data, state)
end


function Gridap.Algebra.solve(
  odeslvr::ODESolver, odeop::ODEOperator,
  t0::Real, tF::Real, us0::Tuple{Vararg{AbstractVector}},
  bhs, yhs, whs, zhs
)
  TSWGenericODESolution(odeslvr, odeop, t0, tF, us0, bhs, yhs, whs, zhs)
end

function Gridap.Algebra.solve(
  odeslvr::ODESolver, odeop::ODEOperator,
  t0::Real, tF::Real, u0::AbstractVector,
  bh, yh, wh, zh
)
  us0 = (u0,)
  bhs = (bh, )
  yhs = (yh, )
  whs = (wh, )
  zhs = (zh, )
  solve(odeslvr, odeop, t0, tF, us0, bhs, yhs, whs, zhs)
end
