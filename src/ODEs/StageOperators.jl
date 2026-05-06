##########################
# DAENonlinearStageOperator #
##########################

struct DAENonlinearStageOperator <: StageOperator
  odeop::ODEOperator
  odeopcache
  tx::Real
  usx::Tuple{Vararg{AbstractVector}}
  ws::Tuple{Vararg{Real}}
end

# NonlinearOperator interface
function Algebra.allocate_residual(
  nlop::DAENonlinearStageOperator, x::AbstractVector
)
  odeop, odeopcache = nlop.odeop, nlop.odeopcache
  tx = nlop.tx
  usx = nlop.usx
  allocate_residual(odeop, tx, usx, odeopcache)
end

function Algebra.residual!(
  r::AbstractVector,
  nlop::DAENonlinearStageOperator, x::AbstractVector
)
  odeop, odeopcache = nlop.odeop, nlop.odeopcache
  tx = nlop.tx
  usx = nlop.usx
  residual!(r, odeop, tx, usx, odeopcache)
end

function Algebra.allocate_jacobian(
  nlop::DAENonlinearStageOperator, x::AbstractVector
)
  odeop, odeopcache = nlop.odeop, nlop.odeopcache
  tx = nlop.tx
  usx = nlop.usx
  allocate_jacobian(odeop, tx, usx, odeopcache)
end

function Algebra.jacobian!(
  J::AbstractMatrix,
  nlop::DAENonlinearStageOperator, x::AbstractVector
)
  odeop, odeopcache = nlop.odeop, nlop.odeopcache
  tx = nlop.tx
  usx = nlop.usx
  ws = nlop.ws
  jacobian!(J, odeop, tx, usx, ws, odeopcache)
  J
end


function Gridap.Algebra.solve!(x::AbstractVector,
  ls::LinearSolver,
  op::DAENonlinearStageOperator,
  cache::Nothing)

  fill!(x,zero(eltype(x)))
  b = residual(op, x)
  A = jacobian(op, x)
  ss = symbolic_setup(ls, A)
  ns = numerical_setup(ss,A)
  rmul!(b,-1)
  solve!(x,ns,b)

  LinearSolverCache(A,b,ns)
end

function Gridap.Algebra.solve!(x::AbstractVector,
  ls::LinearSolver,
  op::DAENonlinearStageOperator,
  cache)

  # println("Non linear stage solve ")
  fill!(x,zero(eltype(x)))
  b = cache.b
  A = cache.A
  ns = cache.ns
  residual!(b, op, x)
  jacobian!(A, op, x)
  numerical_setup!(ns,A)
  rmul!(b,-1)

  solve!(x,ns,b)
  cache
end



function Gridap.Algebra.solve!(x::PVector,
  ls::LinearSolver,
  op::DAENonlinearStageOperator,
  cache::Nothing)

  fill!(x,zero(eltype(x)))
  b = residual(op, x)
  A = jacobian(op, x)
  ss = symbolic_setup(ls, A)
  ns = numerical_setup(ss,A)
  rmul!(b,-1)
  # solve!(x,ns,b)

  _x = Gridap.Algebra.allocate_in_domain(A)
  fill!(_x,0.0)
  Gridap.Algebra.solve!(_x,ns,b)
  copy!(x,_x)
  consistent!(x) |> fetch

  LinearSolverCache(A,b,ns)
end

function Gridap.Algebra.solve!(x::PVector,
  ls::LinearSolver,
  op::DAENonlinearStageOperator,
  cache)

  # println("Non linear stage solve ")
  fill!(x,zero(eltype(x)))
  b = cache.b
  A = cache.A
  ns = cache.ns
  residual!(b, op, x)
  jacobian!(A, op, x)
  numerical_setup!(ns,A)
  rmul!(b,-1)

  # solve!(x,ns,b)
  _x = Gridap.Algebra.allocate_in_domain(A)
  fill!(_x,0.0)
  Gridap.Algebra.solve!(_x,ns,b)
  copy!(x,_x)
  consistent!(x) |> fetch

  cache
end
