
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
  println("my consistent")
  a
end

#### overloads to allow for Pvectors in distributed
function Gridap.Algebra.solve!(
  x::AbstractVector,
  ls::LinearSolver, lop::LinearStageOperator,
  ns::Nothing
)
  println("regular linear stage solver")

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
  println("regular linear stage solver")

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
