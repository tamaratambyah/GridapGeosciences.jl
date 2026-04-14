########################################
# DAEOperator #
########################################
struct TSWOperator <: FEOperator
  progop::FEOperator
  diagop::FEOperator
  diagopb::FEOperator # for consant diagnostics
  diagopw::FEOperator # for bhat
  diagopz::FEOperator # for btilde
  Tfunc::Function # for T = 1/2(h)
  Tspace
  bbarfunc::Function # for ̄bbar = 1/2(bn + bk )
  bbarspace
end

function get_progop(op::TSWOperator)
  op.progop
end

function get_diagop(op::TSWOperator)
  op.diagop
end

function get_diagopb(op::TSWOperator)
  op.diagopb
end

function get_diagopw(op::TSWOperator)
  op.diagopw
end

function get_diagopz(op::TSWOperator)
  op.diagopz
end

function get_Tfunc(op::TSWOperator)
  op.Tfunc
end

function get_Tspace(op::TSWOperator)
  op.Tspace
end

function get_bbarfunc(op::TSWOperator)
  op.bbarfunc
end

function get_bbarspace(op::TSWOperator)
  op.bbarspace
end

function Gridap.ODEs.get_assembler(op::FEOperator)
  op.assem
end

function Gridap.ODEs.get_res(op::FEOperator)
  op.res
end

function Gridap.ODEs.get_jacs(op::FEOperator)
  op.jac
end





########################################
# TransientTSWOperator #
########################################
"""
  const_jac is an optional input relating to the prognostic jacobian.
  Default value is false

  if const_jac = true, the prognostic jacobian is constant over the entire
    simulation, and is computed once at the very beginning of the simulation

  if const_jac = false, the prognostic jacobian is constant over a timestep,
    and is computed at the beginning of each time step
"""
struct TransientTSWOperator <: TransientFEOperator{NonlinearODE}
  op::TSWOperator
  const_jac::Bool
  order::Integer

  function TransientTSWOperator(op::TSWOperator, ;const_jac=false)
    new(op, const_jac, 1)
  end

end


function Gridap.FESpaces.get_algebraic_operator(top::TransientTSWOperator)
  TSWODEOpFromTFEOp(top)
end


function Gridap.FESpaces.get_trial(top::TransientTSWOperator)
  # println("this func")
  progop = get_progop(top.op)
  get_trial(progop)
end



function get_diag_trial(top::TransientTSWOperator)
  # println("diag trial ")
  diagop = get_diagop(top.op)
  get_trial(diagop)
end



function Gridap.Polynomials.get_order(top::TransientTSWOperator)
  top.order
end

function get_progop(top::TransientTSWOperator)
  get_progop(top.op)
end

function get_diagop(top::TransientTSWOperator)
  get_diagop(top.op)
end

function get_diagopb(top::TransientTSWOperator)
  get_diagopb(top.op)
end

function get_diagopw(top::TransientTSWOperator)
  get_diagopw(top.op)
end

function get_diagopz(top::TransientTSWOperator)
  get_diagopz(top.op)
end

function get_Tfunc(top::TransientTSWOperator)
  get_Tfunc(top.op)
end

function get_Tspace(top::TransientTSWOperator)
  get_Tspace(top.op)
end

function get_bbarfunc(top::TransientTSWOperator)
  get_bbarfunc(top.op)
end

function get_bbarspace(top::TransientTSWOperator)
  get_bbarspace(top.op)
end

function is_const_jac(top::TransientTSWOperator)
  top.const_jac
end
