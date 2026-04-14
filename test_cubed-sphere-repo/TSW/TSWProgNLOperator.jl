########################################
# Prognostics  # u, h, B
########################################


struct TSWProgNonlinearOperator <: NonlinearOperator
  odeop
  odeopcache
  ts
  u0
  diagslvr
  J_prog::AbstractMatrix
  const_jac::Bool
  b0
  y
  w
  z
  b
end

# NonlinearOperator interface
function Gridap.Algebra.allocate_residual(
  nlop::TSWProgNonlinearOperator, x::AbstractVector
  )

  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  progop = get_progop(odeop)

  ts = nlop.ts

  diagslvrs = nlop.diagslvr
  diagslvr1,diagslvr2 = diagslvrs


  b0 = nlop.b0
  y = nlop.y
  w = nlop.w
  z = nlop.z
  b = nlop.b

  diagslvrcache = odeopcache.diagslvrcache
  # diagslvrcachew = odeopcache.diagslvrcachew
  diagslvrcachez = odeopcache.diagslvrcachez
  diagslvrcacheb = odeopcache.diagslvrcacheb

  Us = odeopcache.Us
  Ys = odeopcache.Ys
  # Ws = odeopcache.Ws
  Zs = odeopcache.Zs
  Ts = odeopcache.Ts
  Bs = odeopcache.Bs
  Bbars = odeopcache.Bbars

  # solve for b
  diagbop = TSWDiagNonlinearOperatorb(odeop, odeopcache, ts, nlop.u0, x, b0)
  diagslvrcacheb = Gridap.Algebra.solve!(b, diagslvr2, diagbop, diagslvrcacheb)
  update_diagnosticsb!(odeopcache, b, diagslvrcacheb)


  # # Create and solve diagnostics
  diagop = TSWDiagNonlinearOperator(odeop, odeopcache, ts, nlop.u0, x)
  diagslvrcache = Gridap.Algebra.solve!(y, diagslvr1, diagop, diagslvrcache)
  update_diagnostics!(odeopcache, y, diagslvrcache)

  u0 = FEFunction(Us[1], nlop.u0)
  uh = FEFunction(Us[1], x)
  yh = FEFunction(Ys[1], y)
  bh = FEFunction(Bs[1], b)
  bh0 = FEFunction(Bs[1], b0)

  # evaluate T
  Tfunc = get_Tfunc(odeop)
  Th = interpolate(Tfunc(u0,uh), Ts[1])

  # evaluate b bar
  bbarfunc = get_bbarfunc(odeop)
  bbarh = interpolate(bbarfunc(bh0,bh), Bbars[1])

  # # solve for bhat
  # diagopw = TSWDiagNonlinearOperatorw(odeop, odeopcache, ts, b0, b, nlop.u0, x)
  # diagslvrcachew = Gridap.Algebra.solve!(w, diagslvr, diagopw, diagslvrcachew)
  # update_diagnosticsw!(odeopcache, w, diagslvrcachew)
  # wh = FEFunction(Ws[1],w)


  # # solve for btilde
  diagopz = TSWDiagNonlinearOperatorz(odeop, odeopcache, ts, w, b0, b,y)
  diagslvrcachez = Gridap.Algebra.solve!(z, diagslvr1, diagopz, diagslvrcachez)
  update_diagnosticsz!(odeopcache, z, diagslvrcachez)
  zh = FEFunction(Zs[1],z)

  copy_to_cache!(odeopcache, b,y,w,z)

  V = Gridap.FESpaces.get_test(progop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(progop)

  res = Gridap.ODEs.get_res(progop)
  vecdata = collect_cell_vector(V, res(ts, u0, uh, Th, bbarh, zh, yh, v, bh, bh0))
  allocate_vector(assembler, vecdata)


end

function Gridap.Algebra.allocate_jacobian(
  nlop::TSWProgNonlinearOperator, x::AbstractVector
  )
  # println("allocate jacobian")
  # odeop, odeopcache = nlop.odeop, nlop.odeopcache
  # ts = nlop.ts

  # b = odeopcache.diagnosticsb
  # y = odeopcache.diagnostics
  # w = odeopcache.diagnosticsw
  # z = odeopcache.diagnosticsz

  # allocate_prognostic_jacobian(odeop, ts, x, odeopcache, nlop.u0, b, y, w, z)
  nlop.J_prog

end


function Gridap.Algebra.residual!(
  r::AbstractVector,
  nlop::TSWProgNonlinearOperator, x::AbstractVector
  )
  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  progop = get_progop(odeop)


  ts = nlop.ts

  diagslvrs = nlop.diagslvr
  diagslvr1,diagslvr2 = diagslvrs


  b0 = nlop.b0 # odeopcache.diagnostics0 # these have been previously updated
  y = nlop.y # to be solved for
  w = nlop.w # to be solved for
  z = nlop.z # to be solved for
  b = nlop.b

  diagslvrcache = odeopcache.diagslvrcache
  # diagslvrcachew = odeopcache.diagslvrcachew
  diagslvrcachez = odeopcache.diagslvrcachez
  diagslvrcacheb = odeopcache.diagslvrcacheb

  Us = odeopcache.Us
  Ys = odeopcache.Ys
  # Ws = odeopcache.Ws
  Zs = odeopcache.Zs
  Ts = odeopcache.Ts
  Bs = odeopcache.Bs
  Bbars = odeopcache.Bbars


  # solve for b
  diagbop = TSWDiagNonlinearOperatorb(odeop, odeopcache, ts, nlop.u0, x, b0)
  diagslvrcacheb = Gridap.Algebra.solve!(b, diagslvr2, diagbop, diagslvrcacheb)
  update_diagnosticsb!(odeopcache, b, diagslvrcacheb)


  # Create and solve diagnostics
  diagop = TSWDiagNonlinearOperator(odeop, odeopcache, ts, nlop.u0, x)
  diagslvrcache = Gridap.Algebra.solve!(y, diagslvr1, diagop, diagslvrcache)
  update_diagnostics!(odeopcache, y, diagslvrcache)

  u0 = FEFunction(Us[1], nlop.u0)
  uh = FEFunction(Us[1], x)
  yh = FEFunction(Ys[1], y)
  bh0 = FEFunction(Bs[1], b0)
  bh = FEFunction(Bs[1], b)

  # evaluate T
  Tfunc = get_Tfunc(odeop)
  Th = interpolate(Tfunc(u0,uh), Ts[1])

  # evaluate b bar
  bbarfunc = get_bbarfunc(odeop)
  bbarh = interpolate(bbarfunc(bh0,bh), Bbars[1])

  # solve for bhat
  # diagopw = TSWDiagNonlinearOperatorw(odeop, odeopcache, ts, b0, b, nlop.u0, x)
  # diagslvrcachew = Gridap.Algebra.solve!(w, diagslvr, diagopw, diagslvrcachew)
  # update_diagnosticsw!(odeopcache, w, diagslvrcachew)
  # wh = FEFunction(Ws[1],w)

  # solve for btilde
  diagopz = TSWDiagNonlinearOperatorz(odeop, odeopcache, ts, w, b0, b, y)
  diagslvrcachez = Gridap.Algebra.solve!(z, diagslvr1, diagopz, diagslvrcachez)
  update_diagnosticsz!(odeopcache, z, diagslvrcachez)
  zh = FEFunction(Zs[1],z)

  copy_to_cache!(odeopcache, b,y,w,z)

  V = Gridap.FESpaces.get_test(progop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(progop)

  fill!(r, zero(eltype(r)))

  res = Gridap.ODEs.get_res(progop)
  dc = DomainContribution()
  dc = dc + res(ts, u0, uh, Th, bbarh, zh, yh, v, bh, bh0)
  vecdata = collect_cell_vector(V, dc)
  assemble_vector!(r, assembler, vecdata)

  r


end


function Gridap.Algebra.jacobian!(
  J::AbstractMatrix,
  nlop::TSWProgNonlinearOperator, x::AbstractVector
  )


  odeop, odeopcache = nlop.odeop, nlop.odeopcache
  ts = nlop.ts

  b0 = nlop.b0
  b = odeopcache.diagnosticsb
  y = odeopcache.diagnostics
  w = odeopcache.diagnosticsw
  z = odeopcache.diagnosticsz

  Gridap.Algebra.copy_entries!(J, nlop.J_prog)

  if !nlop.const_jac
    prognostic_jacobian!(J, odeop, ts, x, odeopcache, nlop.u0, b, y, w, z, b0)
  end


  J

end



function GridapSolvers.NonlinearSolvers._solve_nr!(x,A,b,dx,ns,nls,op::TSWProgNonlinearOperator)
  # println("using gridap solvers")

  log = nls.log

  # Check for convergence on the initial residual
  res = norm(b)
  done = GridapSolvers.NonlinearSolvers.init!(log,res)

  # Newton-like iterations
  while !done

    # Solve linearized problem
    rmul!(b,-1)
    solve!(dx,ns,b)
    x .+= dx

    # Check convergence for the current residual
    Gridap.Algebra.residual!(b, op, x)
    res  = norm(b)
    done = GridapSolvers.NonlinearSolvers.update!(log,res)

    if !done
      # Update jacobian and solver
      Gridap.Algebra.jacobian!(A, op, x)
      numerical_setup!(ns,A,x)
    end

  end

  GridapSolvers.NonlinearSolvers.finalize!(log,res)
  return x
end
