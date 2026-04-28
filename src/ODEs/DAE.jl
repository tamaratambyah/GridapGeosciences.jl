
########################################
# DAEFEOperator #
########################################

struct DAEFEOperator <: TransientFEOperator{SemilinearODE}
  odefeop::TransientFEOperator
  daefeop::FEOperator
  dae_nl::NonlinearSolver
end

function get_odefeop(op::DAEFEOperator)
  op.odefeop
end

function get_daefeop(op::DAEFEOperator)
  op.daefeop
end

function get_daesolver(op::DAEFEOperator)
  op.dae_nl
end

# TransientFEOperator interface
function Gridap.FESpaces.get_test(op::DAEFEOperator)
  feop = get_odefeop(op)
  Gridap.ODEs.get_test(feop)
end

function Gridap.FESpaces.get_trial(op::DAEFEOperator)
  feop = get_odefeop(op)
  Gridap.ODEs.get_trial(feop)
end

function Gridap.Polynomials.get_order(op::DAEFEOperator)
  feop = get_odefeop(op)
  Gridap.ODEs.get_order(feop)
end

function Gridap.ODEs.get_res(op::DAEFEOperator)
  feop = get_odefeop(op)
  Gridap.ODEs.get_res(feop)
end

function Gridap.ODEs.get_jacs(op::DAEFEOperator)
  feop = get_odefeop(op)
  Gridap.ODEs.get_jacs(feop)
end

function Gridap.ODEs.get_forms(op::DAEFEOperator)
  feop = get_odefeop(op)
  Gridap.ODEs.get_forms(feop)
end

function Gridap.ODEs.is_form_constant(op::DAEFEOperator, k::Integer)
  feop = get_odefeop(op)
  Gridap.ODEs.is_form_constant(feop, k::Integer)
end

function Gridap.ODEs.get_assembler(op::DAEFEOperator)
  feop = get_odefeop(op)
  Gridap.ODEs.get_assembler(feop)
end

### DAEE
function get_res_y(op::DAEFEOperator)
  feop = get_daefeop(op)
  feop.res
end

function get_jacs_y(op::DAEFEOperator)
  feop = get_daefeop(op)
  feop.jac
end

function get_test_y(op::DAEFEOperator)
  feop = get_daefeop(op)
  Gridap.FESpaces.get_test(feop)
end

function get_trial_y(op::DAEFEOperator)
  feop = get_daefeop(op)
  Gridap.FESpaces.get_trial(feop)
end

function get_assembler_y(op::DAEFEOperator)
  feop = get_daefeop(op)
  feop.assem
end

function Gridap.FESpaces.get_algebraic_operator(op::DAEFEOperator)
  DAEODEOpFromTFEOp(op.odefeop,op.daefeop,op.dae_nl)
end

### DAEE FROM WEAK FORM
function get_res_y(op::FEOperatorFromWeakForm)
  op.res
end

function get_jacs_y(op::FEOperatorFromWeakForm)
  op.jac
end

function get_test_y(op::FEOperatorFromWeakForm)
  op.test
end

function get_trial_y(op::FEOperatorFromWeakForm)
  op.trial
end

function get_assembler_y(op::FEOperatorFromWeakForm)
  op.assem
end



######################
# Nonlinear DAE operator #
######################

import Gridap.ODEs: NonlinearOperator
import Gridap.ODEs: ODEOperator
import Gridap.Algebra: LinearSolver
import Gridap.Algebra: LinearSolverCache
import Gridap.Algebra: symbolic_setup, numerical_setup, numerical_setup!
import Gridap.Algebra: NLSolver, NLSolversCache
mutable struct DAENonlinearOperator <: NonlinearOperator
  odeop
  odeopcache
  ui::AbstractVector
  ti::Real
end


function Gridap.Algebra.allocate_residual(
  op::DAENonlinearOperator,x::AbstractVector
)
  zero(x)
end


function Gridap.Algebra.allocate_jacobian(
  disop::DAENonlinearOperator,
  x::AbstractVector
)
  t = disop.ti
  ui = disop.ui

  odeop = disop.odeop
  odeopcache = disop.odeopcache
  Us = odeopcache.Us
  Ys = odeopcache.Ys

  uh = FEFunction(Us[1], ui)

  yh = FEFunction(Ys[1],x)

  Yt = evaluate(get_trial_y(odeop.feop), nothing)
  dy = get_trial_fe_basis(Yt)
  V = get_test_y(odeop.feop)
  v = get_fe_basis(V)
  assembler = get_assembler_y(odeop.feop)

  jac = get_jacs_y(odeop.feop)
  matdata = collect_cell_matrix(Yt, V, jac(t, (uh, yh), dy, v))
  allocate_matrix(assembler, matdata)



end


function Gridap.Algebra.residual!(
  r::AbstractVector,
  disop::DAENonlinearOperator,
  x::AbstractVector
)
  t = disop.ti
  ui = disop.ui

  odeop = disop.odeop
  odeopcache = disop.odeopcache
  Us = odeopcache.Us
  Ys = odeopcache.Ys

  uh = FEFunction(Us[1], ui)

  yh = FEFunction(Ys[1], x)

  V = get_test_y(odeop.feop)
  v = get_fe_basis(V)
  assembler = get_assembler_y(odeop.feop)

  fill!(r, zero(eltype(r)))
  # Residual
  res = get_res_y(odeop.feop)
  vecdata = collect_cell_vector(V, res(t, (uh, yh), v))
  assemble_vector!(r, assembler, vecdata)
  # NLSolver takes a function for J and r, so no need to multiply by -1
  r

end


function Gridap.Algebra.jacobian!(
  J::AbstractMatrix,
  disop::DAENonlinearOperator,
  x::AbstractVector
)
  t = disop.ti
  ui = disop.ui

  odeop = disop.odeop
  odeopcache = disop.odeopcache
  Us = odeopcache.Us
  Ys = odeopcache.Ys

  uh = FEFunction(Us[1], ui)

  yh = FEFunction(Ys[1],x)

  Yt = evaluate(get_trial_y(odeop.feop), nothing)
  dy = get_trial_fe_basis(Yt)
  V = get_test_y(odeop.feop)
  v = get_fe_basis(V)
  assembler = get_assembler_y(odeop.feop)

  jac = get_jacs_y(odeop.feop)
  matdata = collect_cell_matrix(Yt, V, jac(t, (uh, yh), dy, v))
  assemble_matrix!(J, assembler, matdata)
  J


end


function update!(op::DAENonlinearOperator,u::AbstractVector,t::Real)
  copy!(op.ui, u)
  op.ti = t

end

#### Specific solve for DAENonlinearOperator
# want to recompute Jacobian everytime solve! is called, and solve with linear solver
function Gridap.Algebra.solve!(x::AbstractVector,
                ls::LinearSolver,
                op::DAENonlinearOperator,
                cache::Nothing)
  fill!(x,zero(eltype(x)))
  b = Gridap.Algebra.residual(op, x)
  A = Gridap.Algebra.jacobian(op, x)
  ss = symbolic_setup(ls, A)
  ns = numerical_setup(ss,A)
  rmul!(b,-1)

  solve!(x,ns,b)

  LinearSolverCache(A,b,ns)
end

function Gridap.Algebra.solve!(x::AbstractVector,
                ls::LinearSolver,
                op::DAENonlinearOperator,
                cache)
  # println("my DAE linear solver")
  fill!(x,zero(eltype(x)))
  b = cache.b
  A = cache.A
  ns = cache.ns
  Gridap.Algebra.residual!(b, op, x)
  Gridap.Algebra.jacobian!(A, op, x)
  numerical_setup!(ns,A)
  rmul!(b,-1)

  solve!(x,ns,b)

  cache
end

function Gridap.Algebra.solve!(x::PVector,
  ls::LinearSolver,
  op::DAENonlinearOperator,
  cache::Nothing)
  fill!(x,zero(eltype(x)))
  b = Gridap.Algebra.residual(op, x)
  A = Gridap.Algebra.jacobian(op, x)
  ss = symbolic_setup(ls, A)
  ns = numerical_setup(ss,A)
  rmul!(b,-1)

  _x = Gridap.Algebra.allocate_in_domain(A)
  fill!(_x,0.0)
  Gridap.Algebra.solve!(_x,ns,b)
  copy!(x,_x)
  consistent!(x) |> fetch

  LinearSolverCache(A,b,ns)
end

function Gridap.Algebra.solve!(x::PVector,
  ls::LinearSolver,
  op::DAENonlinearOperator,
  cache)
  # println("my DAE linear solver")
  fill!(x,zero(eltype(x)))
  b = cache.b
  A = cache.A
  ns = cache.ns
  Gridap.Algebra.residual!(b, op, x)
  Gridap.Algebra.jacobian!(A, op, x)
  numerical_setup!(ns,A)
  rmul!(b,-1)

  _x = Gridap.Algebra.allocate_in_domain(A)
  fill!(_x,0.0)
  Gridap.Algebra.solve!(_x,ns,b)
  copy!(x,_x)
  consistent!(x) |> fetch

  cache
end

# include a wrapper incase user wants to use a nonlinear solver
function Gridap.Algebra.solve!(x::AbstractVector,nls::NLSolver,op::DAENonlinearOperator,cache::Nothing)
  Gridap.Algebra.solve!(x,nls,op,cache)
end

function Gridap.Algebra.solve!(
  x::AbstractVector,nls::NLSolver,op::DAENonlinearOperator,cache::NLSolversCache)
  # println("my DAE non linear solver")
  Gridap.Algebra.solve!(x,nls,op,cache)
end
