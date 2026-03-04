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
  # solve!(x,ns,b)

  _x = Gridap.Algebra.allocate_in_domain(A)
  fill!(_x,0.0)
  Gridap.Algebra.solve!(_x,ns,b)
  copy!(x,_x)
  consistent!(x) |> fetch

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

  # solve!(x,ns,b)
  _x = Gridap.Algebra.allocate_in_domain(A)
  fill!(_x,0.0)
  Gridap.Algebra.solve!(_x,ns,b)
  copy!(x,_x)
  consistent!(x) |> fetch

  cache
end


#######################
# LinearStageOperator #
#######################


function Gridap.ODEs.LinearStageOperator(
  odeop::DAEODEOperator, odeopcache,
  tx::Real, usx::Tuple{Vararg{AbstractVector}},
  ws::Tuple{Vararg{Real}},
  J::AbstractMatrix, r::AbstractVector, reuse::Bool, sysslvrcache
)
  Gridap.ODEs.residual!(r, odeop, tx, usx, odeopcache)

  if isnothing(sysslvrcache) || !reuse
    Gridap.ODEs.jacobian!(J, odeop, tx, usx, ws, odeopcache)
  end

  LinearStageOperator(J, r, tx, ws, usx, reuse)
end

## helper for nonPvectors
function PartitionedArrays.consistent!(a::Vector)
  # println("my consistent")
  a
end

#### overloads to allow for Pvectors in distributed
function Gridap.Algebra.solve!(
  x::AbstractVector,
  ls::LinearSolver, lop::LinearStageOperator,
  ns::Nothing
)
  # println("regular linear stage solver")

  J = lop.J
  ss = symbolic_setup(ls, J)
  ns = numerical_setup(ss, J)

  r = lop.r
  rmul!(r, -1)

  # solve!(x, ns, r)
  _x = Gridap.Algebra.allocate_in_domain(J)
  fill!(_x,0.0)
  Gridap.Algebra.solve!(_x,ns,r)
  copy!(x,_x)
  consistent!(x) |> fetch

  ns
end

function Gridap.Algebra.solve!(
  x::AbstractVector,
  ls::LinearSolver, lop::LinearStageOperator,
  ns
)
  # println("regular linear stage solver")

  J = lop.J
  if !lop.reuse
    J = lop.J
    numerical_setup!(ns, J)
  end

  r = lop.r
  rmul!(r, -1)

  # solve!(x, ns, r)
  _x = Gridap.Algebra.allocate_in_domain(J)
  fill!(_x,0.0)
  Gridap.Algebra.solve!(_x,ns,r)
  copy!(x,_x)
  consistent!(x) |> fetch

  ns
end
