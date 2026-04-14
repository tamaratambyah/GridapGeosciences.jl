

########################################
# Diagnostics  # F, Φ, q
########################################

struct TSWDiagNonlinearOperator <: NonlinearOperator
  odeop
  odeopcache
  ts
  u0
  u
end


# NonlinearOperator interface
function Gridap.Algebra.allocate_residual(
  nlop::TSWDiagNonlinearOperator, x::AbstractVector
  )

  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagop = get_diagop(odeop)

  ts = nlop.ts

  Us = odeopcache.Us
  Ys = odeopcache.Ys

  u0 = EvaluationFunction(Us[1], nlop.u0)
  uh = EvaluationFunction(Us[1], nlop.u)
  yh = EvaluationFunction(Ys[1], x)

  V = Gridap.FESpaces.get_test(diagop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagop)

  res = Gridap.ODEs.get_res(diagop)
  vecdata = collect_cell_vector(V, res(ts, u0, uh,  yh, v))
  allocate_vector(assembler, vecdata)



end

function Gridap.Algebra.allocate_jacobian(
  nlop::TSWDiagNonlinearOperator, x::AbstractVector
  )
  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagop = get_diagop(odeop)

  ts = nlop.ts

  Us = odeopcache.Us
  Ys = odeopcache.Ys

  u0 = EvaluationFunction(Us[1], nlop.u0)
  uh = EvaluationFunction(Us[1], nlop.u)

  yh = EvaluationFunction(Ys[1], x)

  Yt = evaluate(Gridap.FESpaces.get_trial(diagop), nothing)
  dy = get_trial_fe_basis(Yt)
  V = Gridap.FESpaces.get_test(diagop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagop)

  jac = Gridap.ODEs.get_jacs(diagop)
  dc = DomainContribution()
  dc = dc + jac(ts, u0, uh, yh, dy, v)

  matdata = collect_cell_matrix(Yt, V, dc)
  allocate_matrix(assembler, matdata)


end


function Gridap.Algebra.residual!(
  r::AbstractVector,
  nlop::TSWDiagNonlinearOperator, x::AbstractVector
  )
  # println("disagnostic solver")
  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagop = get_diagop(odeop)

  ts = nlop.ts

  Us = odeopcache.Us
  Ys = odeopcache.Ys

  u0 = EvaluationFunction(Us[1], nlop.u0)
  uh = EvaluationFunction(Us[1], nlop.u)
  yh = EvaluationFunction(Ys[1], x)

  V = Gridap.FESpaces.get_test(diagop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagop)

  fill!(r, zero(eltype(r)))

  res = Gridap.ODEs.get_res(diagop)
  dc = DomainContribution()
  dc = dc + res(ts, u0, uh, yh, v)
  vecdata = collect_cell_vector(V, dc)
  assemble_vector!(r, assembler, vecdata)

  r

end



function Gridap.Algebra.jacobian!(
  J::AbstractMatrix,
  nlop::TSWDiagNonlinearOperator, x::AbstractVector
  )


  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagop = get_diagop(odeop)

  ts = nlop.ts

  Us = odeopcache.Us
  Ys = odeopcache.Ys

  u0 = EvaluationFunction(Us[1], nlop.u0)
  uh = EvaluationFunction(Us[1], nlop.u)

  yh = EvaluationFunction(Ys[1], x)

  Yt = evaluate(Gridap.FESpaces.get_trial(diagop), nothing)
  dy = get_trial_fe_basis(Yt)
  V = Gridap.FESpaces.get_test(diagop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagop)

  jac = Gridap.ODEs.get_jacs(diagop)
  dc = DomainContribution()
  dc = dc + jac(ts, u0, uh, yh, dy, v)

  matdata = collect_cell_matrix(Yt, V, dc)
  assemble_matrix!(J, assembler, matdata)

  J

end



#### Specific solve for TSWDiagNonlinearOperator
# want to recompute Jacobian everytime solve! is called, and solve with linear solver
function Gridap.Algebra.solve!(x::AbstractVector,
                ls::LinearSolver,
                op::TSWDiagNonlinearOperator,
                cache::Nothing)
  # println("diagnostic solve 0")

  fill!(x,zero(eltype(x)))
  b = Gridap.Algebra.residual(op, x)
  A = Gridap.Algebra.jacobian(op, x)
  ss = symbolic_setup(ls, A)
  ns = numerical_setup(ss,A)
  rmul!(b,-1)
  Gridap.Algebra.solve!(x,ns,b)
  LinearSolverCache(A,b,ns)
end

function Gridap.Algebra.solve!(x::AbstractVector,
                ls::LinearSolver,
                op::TSWDiagNonlinearOperator,
                cache)

  # println("diagnostic solve")

  fill!(x,zero(eltype(x)))

  b = cache.b
  A = cache.A
  ns = cache.ns
  Gridap.Algebra.residual!(b, op, x)
  Gridap.Algebra.jacobian!(A, op, x)
  numerical_setup!(ns,A)
  rmul!(b,-1)
  Gridap.Algebra.solve!(x,ns,b)
  cache
end






########################################
# Initial diagnostics  # b
# solve once per time step
########################################

struct TSWDiagNonlinearOperatorb <: NonlinearOperator
  odeop
  odeopcache
  ts
  u0
  u
  b0
end

# NonlinearOperator interface
function Gridap.Algebra.allocate_residual(
  nlop::TSWDiagNonlinearOperatorb, x::AbstractVector
)

  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagopb = get_diagopb(odeop)

  ts = nlop.ts

  Us = odeopcache.Us
  Bs = odeopcache.Bs

  u0 = EvaluationFunction(Us[1], nlop.u0)
  uh = EvaluationFunction(Us[1], nlop.u)
  bh = EvaluationFunction(Bs[1], x)
  bh0 = EvaluationFunction(Bs[1], nlop.b0)


  V = Gridap.FESpaces.get_test(diagopb)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagopb)

  res = Gridap.ODEs.get_res(diagopb)
  vecdata = collect_cell_vector(V, res(ts, u0, uh, bh, v, bh0))
  allocate_vector(assembler, vecdata)

end

function Gridap.Algebra.allocate_jacobian(
  nlop::TSWDiagNonlinearOperatorb, x::AbstractVector
)
  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagopb = get_diagopb(odeop)

  ts = nlop.ts

  Us = odeopcache.Us
  Bs = odeopcache.Bs

  u0 = EvaluationFunction(Us[1], nlop.u0)
  uh = EvaluationFunction(Us[1], nlop.u)
  bh = EvaluationFunction(Bs[1], x)
  bh0 = EvaluationFunction(Bs[1], nlop.b0)

  Yt0 = evaluate(Gridap.FESpaces.get_trial(diagopb), nothing)
  db = get_trial_fe_basis(Yt0)
  V = Gridap.FESpaces.get_test(diagopb)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagopb)

  jac = Gridap.ODEs.get_jacs(diagopb)
  dc = DomainContribution()
  dc = dc + jac(ts, u0, uh, bh, db, v, bh0)

  matdata = collect_cell_matrix(Yt0, V, dc)
  allocate_matrix(assembler, matdata)


end


function Gridap.Algebra.residual!(
  r::AbstractVector,
  nlop::TSWDiagNonlinearOperatorb, x::AbstractVector
)
  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagopb = get_diagopb(odeop)

  ts = nlop.ts

  Us = odeopcache.Us
  Bs = odeopcache.Bs

  u0 = EvaluationFunction(Us[1], nlop.u0)
  uh = EvaluationFunction(Us[1], nlop.u)
  bh = EvaluationFunction(Bs[1], x)
  bh0 = EvaluationFunction(Bs[1], nlop.b0)


  V = Gridap.FESpaces.get_test(diagopb)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagopb)

  fill!(r, zero(eltype(r)))

  res = Gridap.ODEs.get_res(diagopb)
  dc = DomainContribution()
  dc = dc + res(ts, u0, uh, bh, v, bh0)
  vecdata = collect_cell_vector(V, dc)
  assemble_vector!(r, assembler, vecdata)

  r

end



function Gridap.Algebra.jacobian!(
  J::AbstractMatrix,
  nlop::TSWDiagNonlinearOperatorb, x::AbstractVector
)


  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagopb = get_diagopb(odeop)

  ts = nlop.ts

  Us = odeopcache.Us
  Bs = odeopcache.Bs

  u0 = EvaluationFunction(Us[1], nlop.u0)
  uh = EvaluationFunction(Us[1], nlop.u)
  bh = EvaluationFunction(Bs[1], x)
  bh0 = EvaluationFunction(Bs[1], nlop.b0)


  Yt0 = evaluate(Gridap.FESpaces.get_trial(diagopb), nothing)
  db = get_trial_fe_basis(Yt0)
  V = Gridap.FESpaces.get_test(diagopb)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagopb)

  jac = Gridap.ODEs.get_jacs(diagopb)
  dc = DomainContribution()
  dc = dc + jac(ts, u0, uh, bh, db, v, bh0)

  matdata = collect_cell_matrix(Yt0, V, dc)
  assemble_matrix!(J, assembler, matdata)

  J

end


### Specific solve for TSWDiagNonlinearOperatorb
# want to recompute Jacobian everytime solve! is called, and solve with linear solver
function Gridap.Algebra.solve!(x::AbstractVector,
  ls::LinearSolver,
  op::TSWDiagNonlinearOperatorb,
  cache::Nothing)
# println("diagnostic solve 0")

  fill!(x,zero(eltype(x)))
  b = Gridap.Algebra.residual(op, x)
  A = Gridap.Algebra.jacobian(op, x)
  ss = symbolic_setup(ls, A)
  ns = numerical_setup(ss,A)
  rmul!(b,-1)
  Gridap.Algebra.solve!(x,ns,b)
  LinearSolverCache(A,b,ns)
end

function Gridap.Algebra.solve!(x::AbstractVector,
  ls::LinearSolver,
  op::TSWDiagNonlinearOperatorb,
  cache)

  # println("diagnostic solve")

  fill!(x,zero(eltype(x)))

  b = cache.b
  A = cache.A
  ns = cache.ns
  Gridap.Algebra.residual!(b, op, x)
  Gridap.Algebra.jacobian!(A, op, x)
  numerical_setup!(ns,A)
  rmul!(b,-1)
  Gridap.Algebra.solve!(x,ns,b)
  cache
end





########################################
# Diagnostics for bhat # b2
########################################

struct TSWDiagNonlinearOperatorw <: NonlinearOperator
  odeop
  odeopcache
  ts
  b0
  b
  u0
  u
end

# NonlinearOperator interface
function Gridap.Algebra.allocate_residual(
  nlop::TSWDiagNonlinearOperatorw, x::AbstractVector
)

  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagopw = get_diagopw(odeop)

  ts = nlop.ts

  Bs = odeopcache.Bs
  Ws = odeopcache.Ws
  Us = odeopcache.Us

  u0 = EvaluationFunction(Us[1],nlop.u0)
  uh = EvaluationFunction(Us[1],nlop.u)
  bh0 = EvaluationFunction(Bs[1], nlop.b0)
  bh = EvaluationFunction(Bs[1], nlop.b)
  wh = EvaluationFunction(Ws[1], x)

  V = Gridap.FESpaces.get_test(diagopw)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagopw)

  res = Gridap.ODEs.get_res(diagopw)
  vecdata = collect_cell_vector(V, res(ts, u0,uh,bh0, bh, wh, v))
  allocate_vector(assembler, vecdata)

end

function Gridap.Algebra.allocate_jacobian(
  nlop::TSWDiagNonlinearOperatorw, x::AbstractVector
)
  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagopw = get_diagopw(odeop)

  ts = nlop.ts

  Ws = odeopcache.Ws
  Us = odeopcache.Us

  u0 = EvaluationFunction(Us[1],nlop.u0)
  uh = EvaluationFunction(Us[1],nlop.u)
  wh = EvaluationFunction(Ws[1], x)


  Wt = evaluate(Gridap.FESpaces.get_trial(diagopw), nothing)
  dw = get_trial_fe_basis(Wt)
  V = Gridap.FESpaces.get_test(diagopw)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagopw)

  jac = Gridap.ODEs.get_jacs(diagopw)
  dc = DomainContribution()
  dc = dc + jac(ts, u0, uh, wh, dw, v)

  matdata = collect_cell_matrix(Wt, V, dc)
  allocate_matrix(assembler, matdata)


end


function Gridap.Algebra.residual!(
  r::AbstractVector,
  nlop::TSWDiagNonlinearOperatorw, x::AbstractVector
)
  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagopw = get_diagopw(odeop)

  ts = nlop.ts

  Bs = odeopcache.Bs
  Ws = odeopcache.Ws
  Us = odeopcache.Us

  u0 = EvaluationFunction(Us[1],nlop.u0)
  uh = EvaluationFunction(Us[1],nlop.u)
  bh0 = EvaluationFunction(Bs[1], nlop.b0)
  bh = EvaluationFunction(Bs[1], nlop.b)
  wh = EvaluationFunction(Ws[1], x)


  V = Gridap.FESpaces.get_test(diagopw)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagopw)

  fill!(r, zero(eltype(r)))

  res = Gridap.ODEs.get_res(diagopw)
  dc = DomainContribution()
  dc = dc + res(ts, u0, uh, bh0, bh, wh, v)
  vecdata = collect_cell_vector(V, dc)
  assemble_vector!(r, assembler, vecdata)

  r

end



function Gridap.Algebra.jacobian!(
  J::AbstractMatrix,
  nlop::TSWDiagNonlinearOperatorw, x::AbstractVector
)


  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagopw = get_diagopw(odeop)

  ts = nlop.ts

  Ws = odeopcache.Ws
  Us = odeopcache.Us

  u0 = EvaluationFunction(Us[1],nlop.u0)
  uh = EvaluationFunction(Us[1],nlop.u)
  wh = EvaluationFunction(Ws[1], x)


  Wt = evaluate(Gridap.FESpaces.get_trial(diagopw), nothing)
  dw = get_trial_fe_basis(Wt)
  V = Gridap.FESpaces.get_test(diagopw)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagopw)

  jac = Gridap.ODEs.get_jacs(diagopw)
  dc = DomainContribution()
  dc = dc + jac(ts, u0, uh, wh, dw, v)

  matdata = collect_cell_matrix(Wt, V, dc)
  assemble_matrix!(J, assembler, matdata)

  J

end

## LHS matrix is not weighted, so do not specific solve. Just use general one
#### Specific solve for TSWDiagNonlinearOperatorz
# want to recompute Jacobian everytime solve! is called, and solve with linear solver
function Gridap.Algebra.solve!(x::AbstractVector,
  ls::LinearSolver,
  op::TSWDiagNonlinearOperatorw,
  cache::Nothing)
# println("diagnostic solve 0")

  fill!(x,zero(eltype(x)))
  b = Gridap.Algebra.residual(op, x)
  A = Gridap.Algebra.jacobian(op, x)
  ss = symbolic_setup(ls, A)
  ns = numerical_setup(ss,A)
  rmul!(b,-1)
  Gridap.Algebra.solve!(x,ns,b)
  LinearSolverCache(A,b,ns)
end

function Gridap.Algebra.solve!(x::AbstractVector,
  ls::LinearSolver,
  op::TSWDiagNonlinearOperatorw,
  cache)

  # println("diagnostic solve")

  fill!(x,zero(eltype(x)))

  b = cache.b
  A = cache.A
  ns = cache.ns
  Gridap.Algebra.residual!(b, op, x)
  Gridap.Algebra.jacobian!(A, op, x)
  numerical_setup!(ns,A)
  rmul!(b,-1)
  Gridap.Algebra.solve!(x,ns,b)
  cache
end

########################################
# Diagnostics for  btilde # b3
########################################



struct TSWDiagNonlinearOperatorz <: NonlinearOperator
  odeop
  odeopcache
  ts
  w
  b0
  b
  y
end

# NonlinearOperator interface
function Gridap.Algebra.allocate_residual(
  nlop::TSWDiagNonlinearOperatorz, x::AbstractVector
)

  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagopz = get_diagopz(odeop)

  ts = nlop.ts


  # Ws = odeopcache.Ws
  Zs = odeopcache.Zs
  Bs = odeopcache.Bs
  Bbars = odeopcache.Bbars
  Ys = odeopcache.Ys

  bh = EvaluationFunction(Bs[1], nlop.b)
  bh0 = EvaluationFunction(Bs[1], nlop.b0)
  # wh = EvaluationFunction(Ws[1], nlop.w)
  zh = EvaluationFunction(Zs[1], x)
  yh = EvaluationFunction(Ys[1], nlop.y)

  # evaluate b bar
  bbarfunc = get_bbarfunc(odeop)
  bbarh = interpolate(bbarfunc(bh0,bh), Bbars[1])

  V = Gridap.FESpaces.get_test(diagopz)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagopz)

  res = Gridap.ODEs.get_res(diagopz)
  vecdata = collect_cell_vector(V, res(ts, yh, bbarh, zh, v, bh, bh0))
  allocate_vector(assembler, vecdata)

end

function Gridap.Algebra.allocate_jacobian(
  nlop::TSWDiagNonlinearOperatorz, x::AbstractVector
)
  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagopz = get_diagopz(odeop)

  ts = nlop.ts


  # Ws = odeopcache.Ws
  Zs = odeopcache.Zs
  Bs = odeopcache.Bs
  Bbars = odeopcache.Bbars
  Ys = odeopcache.Ys

  bh = EvaluationFunction(Bs[1], nlop.b)
  bh0 = EvaluationFunction(Bs[1], nlop.b0)
  # wh = EvaluationFunction(Ws[1], nlop.w)
  zh = EvaluationFunction(Zs[1], x)
  yh = EvaluationFunction(Ys[1], nlop.y)

  # evaluate b bar
  bbarfunc = get_bbarfunc(odeop)
  bbarh = interpolate(bbarfunc(bh0,bh), Bbars[1])

  Zt = evaluate(Gridap.FESpaces.get_trial(diagopz), nothing)
  dz = get_trial_fe_basis(Zt)
  V = Gridap.FESpaces.get_test(diagopz)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagopz)

  jac = Gridap.ODEs.get_jacs(diagopz)
  dc = DomainContribution()
  dc = dc + jac(ts, yh,bbarh, zh, dz, v, bh)

  matdata = collect_cell_matrix(Zt, V, dc)
  allocate_matrix(assembler, matdata)


end


function Gridap.Algebra.residual!(
  r::AbstractVector,
  nlop::TSWDiagNonlinearOperatorz, x::AbstractVector
)
  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagopz = get_diagopz(odeop)

  ts = nlop.ts


  # Ws = odeopcache.Ws
  Zs = odeopcache.Zs
  Bs = odeopcache.Bs
  Bbars = odeopcache.Bbars
  Ys = odeopcache.Ys

  bh = EvaluationFunction(Bs[1], nlop.b)
  bh0 = EvaluationFunction(Bs[1], nlop.b0)
  # wh = EvaluationFunction(Ws[1], nlop.w)
  zh = EvaluationFunction(Zs[1], x)
  yh = EvaluationFunction(Ys[1], nlop.y)

  # evaluate b bar
  bbarfunc = get_bbarfunc(odeop)
  bbarh = interpolate(bbarfunc(bh0,bh), Bbars[1])


  V = Gridap.FESpaces.get_test(diagopz)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagopz)

  fill!(r, zero(eltype(r)))

  res = Gridap.ODEs.get_res(diagopz)
  dc = DomainContribution()
  dc = dc + res(ts, yh, bbarh, zh, v,bh, bh0)
  vecdata = collect_cell_vector(V, dc)
  assemble_vector!(r, assembler, vecdata)

  r

end



function Gridap.Algebra.jacobian!(
  J::AbstractMatrix,
  nlop::TSWDiagNonlinearOperatorz, x::AbstractVector
)


  odeop, odeopcache = nlop.odeop, nlop.odeopcache

  diagopz = get_diagopz(odeop)

  ts = nlop.ts


  # Ws = odeopcache.Ws
  Zs = odeopcache.Zs
  Bs = odeopcache.Bs
  Bbars = odeopcache.Bbars
  Ys = odeopcache.Ys

  bh = EvaluationFunction(Bs[1], nlop.b)
  bh0 = EvaluationFunction(Bs[1], nlop.b0)
  # wh = EvaluationFunction(Ws[1], nlop.w)
  zh = EvaluationFunction(Zs[1], x)
  yh = EvaluationFunction(Ys[1], nlop.y)

  # evaluate b bar
  bbarfunc = get_bbarfunc(odeop)
  bbarh = interpolate(bbarfunc(bh0,bh), Bbars[1])


  Zt = evaluate(Gridap.FESpaces.get_trial(diagopz), nothing)
  dz = get_trial_fe_basis(Zt)
  V = Gridap.FESpaces.get_test(diagopz)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(diagopz)

  jac = Gridap.ODEs.get_jacs(diagopz)
  dc = DomainContribution()
  dc = dc + jac(ts, yh, bbarh, zh, dz, v,bh)

  matdata = collect_cell_matrix(Zt, V, dc)
  assemble_matrix!(J, assembler, matdata)

  J

end


#### Specific solve for TSWDiagNonlinearOperatorz
# want to recompute Jacobian everytime solve! is called, and solve with linear solver
function Gridap.Algebra.solve!(x::AbstractVector,
  ls::LinearSolver,
  op::TSWDiagNonlinearOperatorz,
  cache::Nothing)
# println("diagnostic solve 0")

  fill!(x,zero(eltype(x)))
  b = Gridap.Algebra.residual(op, x)
  A = Gridap.Algebra.jacobian(op, x)
  ss = symbolic_setup(ls, A)
  ns = numerical_setup(ss,A)
  rmul!(b,-1)
  Gridap.Algebra.solve!(x,ns,b)
  LinearSolverCache(A,b,ns)
end

function Gridap.Algebra.solve!(x::AbstractVector,
  ls::LinearSolver,
  op::TSWDiagNonlinearOperatorz,
  cache)

  # println("diagnostic solve")

  fill!(x,zero(eltype(x)))

  b = cache.b
  A = cache.A
  ns = cache.ns
  Gridap.Algebra.residual!(b, op, x)
  Gridap.Algebra.jacobian!(A, op, x)
  numerical_setup!(ns,A)
  rmul!(b,-1)
  Gridap.Algebra.solve!(x,ns,b)
  cache
end
