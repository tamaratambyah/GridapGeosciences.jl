import Gridap.ODEs: get_order
import Gridap.ODEs: get_trial
import Gridap.ODEs: get_test
import Gridap.ODEs: allocate_space
import Gridap.ODEs: _make_uh_from_us
import Gridap.ODEs: get_jacs
import Gridap.FESpaces: collect_cell_vector
import Gridap.FESpaces: collect_cell_matrix
import Gridap.FESpaces: assemble_vector_add!
import Gridap.FESpaces: allocate_matrix
import Gridap.FESpaces: assemble_matrix_add!
import Gridap.ODEs: AbstractLinearODE
import Gridap.ODEs: AbstractQuasilinearODE
import Gridap.ODEs: AbstractSemilinearODE
import Gridap.Algebra: allocate_vector
import Gridap.ODEs: ODEOperatorType
import Gridap.Helpers: @unreachable
import Gridap.ODEs: allocate_tfeopcache
import Gridap.CellData: DomainContribution, num_domains

#######################
# DAEODEOpFromTFEOpCache #
#######################
"""
    struct DAEODEOpFromTFEOpCache <: GridapType

Structure that stores the `TransientFESpace` and cache of a
`TransientFEOperator`, as well as the jacobian matrices and residual if they
are constant.
"""
mutable struct DAEODEOpFromTFEOpCache <: GridapType
  Us
  Uts
  tfeopcache
  const_forms
  diagnostics
  diagnostics0
  Ys
  Yts
end

##################
# ODEOpFromTFEOp #
##################
"""
    struct ODEOpFromTFEOp <: ODEOperator end

Wrapper that transforms a `TransientFEOperator` into an `ODEOperator`, i.e.
takes `residual(t, uh, ∂t[uh], ..., ∂t^N[uh], vh)` and returns
`residual(t, us)`, where `us[k] = ∂t^k[us]` and `uf` represents the free values
of `uh`.
"""
struct DAEODEOpFromTFEOp{T} <: DAEODEOperator{T}
  tfeop::TransientFEOperator{T}
  feop::FEOperator
  dae_nl::NonlinearSolver

  function DAEODEOpFromTFEOp(tfeop::TransientFEOperator{T},feop::FEOperator,dae_nl::NonlinearSolver) where {T}
    order = get_order(tfeop)
    if order == 0
      is_quasilinear = T <: AbstractQuasilinearODE
      is_linear = T <: AbstractLinearODE
      if is_quasilinear && !is_linear
        msg = """
        For an operator of order zero, the definitions of quasilinear,
        semilinear and linear coincide. Make sure that you have defined the
        transient FE operator as linear.
        """
        @unreachable msg
      else
        new{T}(tfeop, feop, dae_nl)
      end
    else
      new{T}(tfeop, feop, dae_nl)
    end
  end
end

# ODEOperator interface
function Gridap.Polynomials.get_order(odeop::DAEODEOpFromTFEOp)
  Gridap.ODEs.get_order(odeop.tfeop)
end

function Gridap.ODEs.get_num_forms(odeop::DAEODEOpFromTFEOp)
  Gridap.ODEs.get_num_forms(odeop.tfeop)
end

function Gridap.ODEs.get_forms(odeop::DAEODEOpFromTFEOp)
  Gridap.ODEs.get_forms(odeop.tfeop)
end

function Gridap.ODEs.is_form_constant(odeop::DAEODEOpFromTFEOp, k::Integer)
  Gridap.ODEs.is_form_constant(odeop.tfeop, k)
end

function Gridap.ODEs.allocate_odeopcache(
  odeop::DAEODEOpFromTFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}}
)
  println("allocating my cache")

  ### DIAGNOSTICS

  # Allocate FE spaces for derivatives - diagnostics
  order = get_order(odeop)
  Yt = get_trial_y(odeop.feop)
  Y = allocate_space(Yt)
  Yts = (Yt,)
  Ys = (Y,)
  for k in 1:order
    Yts = (Yts..., ∂t(Yts[k]))
    Ys = (Ys..., allocate_space(Yts[k+1]))
  end

  # set initial diagnostics to be uh
  diagnostics = zero(Ys[1])
  diagnostics0 = zero(Ys[1])
  println("set diagnostics")

  # Allocate FE spaces for derivatives
  order = get_order(odeop)
  Ut = get_trial(odeop.tfeop)
  U = allocate_space(Ut)
  Uts = (Ut,)
  Us = (U,)
  for k in 1:order
    Uts = (Uts..., ∂t(Uts[k]))
    Us = (Us..., allocate_space(Uts[k+1]))
  end

  # Allocate the cache of the FE operator
  tfeopcache = allocate_tfeopcache(odeop.tfeop, t, us)

  # Variables for assembly
  uh = _make_uh_from_us(odeop, us, Us)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  Ut = evaluate(get_trial(odeop.tfeop), nothing)
  du = get_trial_fe_basis(Ut)
  assembler = Gridap.ODEs.get_assembler(odeop.tfeop)

  # Store the forms that are constant
  const_forms = ()
  num_forms = Gridap.ODEs.get_num_forms(odeop.tfeop)
  jacs = get_jacs(odeop.tfeop)

  yh = zero(Ys[1])
  yh0 = zero(Ys[1])
  # We want the stored jacobians to have the same sparsity as the full jacobian
  # (when all orders are considered), so we start by allocating it and we will assemble
  # the constant jacobians in a copy of the full jacobian
  # We need a little workaround here since when the `ODEOperator` is quasilinear or
  # semilinear but not linear, it has only one form but `order+1` jacobians.
  dc = DomainContribution()
  for k in 0:order
    jac = jacs[k+1]
    dc = dc + jac(t, (uh,yh), du, v,yh0)
  end
  matdata = collect_cell_matrix(Ut, V, dc)
  J_full = allocate_matrix(assembler, matdata)

  odeoptype = ODEOperatorType(odeop)
  if odeoptype <: AbstractLinearODE
    for k in 0:num_forms-1
      const_form = nothing
      if Gridap.ODEs.is_form_constant(odeop, k)
        jac = jacs[k+1]
        dc = jac(t, (uh,yh), du, v,yh0)
        matdata = collect_cell_matrix(Ut, V, dc)
        const_form = copy(J_full)
        LinearAlgebra.fillstored!(const_form, zero(eltype(const_form)))
        assemble_matrix_add!(const_form, assembler, matdata)
      end
      const_forms = (const_forms..., const_form)
    end
  elseif odeoptype <: AbstractQuasilinearODE
    const_form = nothing
    k = order
    if Gridap.ODEs.is_form_constant(odeop, k)
      jac = jacs[k+1]
      dc = jac(t, (uh,yh), du, v,yh0)
      matdata = collect_cell_matrix(Ut, V, dc)
      const_form = copy(J_full)
      LinearAlgebra.fillstored!(const_form, zero(eltype(const_form)))
      assemble_matrix_add!(const_form, assembler, matdata)
    end
    const_forms = (const_forms..., const_form)
  end



  DAEODEOpFromTFEOpCache(Us, Uts, tfeopcache, const_forms,diagnostics,diagnostics0,Ys,Yts)
end

function Gridap.ODEs.update_odeopcache!(odeopcache, odeop::DAEODEOpFromTFEOp, t::Real,diagnostics,diagnostics0)
  println("updating my cache with initial diagnostics")

  Us = ()
  for k in 0:get_order(odeop)
    Us = (Us..., evaluate!(odeopcache.Us[k+1], odeopcache.Uts[k+1], t))
  end
  odeopcache.Us = Us

  tfeopcache, tfeop = odeopcache.tfeopcache, odeop.tfeop
  odeopcache.tfeopcache = update_tfeopcache!(tfeopcache, tfeop, t)

  # update diagnostics
  Ys = ()
  for k in 0:get_order(odeop)
    Ys = (Ys..., evaluate!(odeopcache.Ys[k+1], odeopcache.Yts[k+1], t))
  end
  odeopcache.Ys = Ys

  odeopcache.diagnostics = FEFunction(Ys[1],diagnostics)
  odeopcache.diagnostics0 = FEFunction(Ys[1],diagnostics0)

  odeopcache
end

function Gridap.ODEs.update_odeopcache!(odeopcache, odeop::DAEODEOpFromTFEOp, t::Real,diagnostics)
  println("updating my cache")

  Us = ()
  for k in 0:get_order(odeop)
    Us = (Us..., evaluate!(odeopcache.Us[k+1], odeopcache.Uts[k+1], t))
  end
  odeopcache.Us = Us

  tfeopcache, tfeop = odeopcache.tfeopcache, odeop.tfeop
  odeopcache.tfeopcache = update_tfeopcache!(tfeopcache, tfeop, t)

  # update diagnostics
  Ys = ()
  for k in 0:get_order(odeop)
    Ys = (Ys..., evaluate!(odeopcache.Ys[k+1], odeopcache.Yts[k+1], t))
  end
  odeopcache.Ys = Ys

  odeopcache.diagnostics = FEFunction(Ys[1],diagnostics)

  odeopcache
end

#####
function Gridap.Algebra.allocate_residual(
  odeop::DAEODEOpFromTFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(odeop.tfeop)

  yh = odeopcache.diagnostics
  yh0 = odeopcache.diagnostics0

  res = Gridap.ODEs.get_res(odeop.tfeop)
  vecdata = collect_cell_vector(V, res(t, (uh,yh), v,yh0))
  allocate_vector(assembler, vecdata)
end

function Gridap.Algebra.residual!(
  r::AbstractVector, odeop::DAEODEOpFromTFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache; add::Bool=false
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(odeop.tfeop)

  !add && fill!(r, zero(eltype(r)))

  yh = odeopcache.diagnostics
  yh0 = odeopcache.diagnostics0

  res = Gridap.ODEs.get_res(odeop.tfeop)
  dc = res(t, (uh,yh), v, yh0)
  vecdata = collect_cell_vector(V, dc)
  assemble_vector_add!(r, assembler, vecdata)

  r
end

function Gridap.Algebra.residual!(
  r::AbstractVector, odeop::DAEODEOpFromTFEOp{<:AbstractQuasilinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache; add::Bool=false
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(odeop.tfeop)

  !add && fill!(r, zero(eltype(r)))

  yh = odeopcache.diagnostics
  yh0 = odeopcache.diagnostics0

  # Residual
  res = Gridap.ODEs.get_res(odeop.tfeop)
  dc = res(t, (uh,yh), v, yh0)

  # Mass
  order = get_order(odeop)
  mass = Gridap.ODEs.get_forms(odeop.tfeop)[1]
  ∂tNuh = ∂t(uh, Val(order))
  dc = dc + mass(t, uh, ∂tNuh, v)

  vecdata = collect_cell_vector(V, dc)
  assemble_vector_add!(r, assembler, vecdata)

  r
end

function Gridap.Algebra.residual!(
  r::AbstractVector, odeop::DAEODEOpFromTFEOp{<:AbstractSemilinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache; add::Bool=false
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(odeop.tfeop)

  !add && fill!(r, zero(eltype(r)))

  yh = odeopcache.diagnostics
  yh0 = odeopcache.diagnostics0

  # Residual
  res = Gridap.ODEs.get_res(odeop.tfeop)
  dc = res(t, (uh,yh), v, yh0)

  # Mass
  order = get_order(odeop)
  mass = Gridap.ODEs.get_forms(odeop.tfeop)[1]
  ∂tNuh = ∂t(uh, Val(order))
  dc = dc + mass(t, ∂tNuh, v)

  vecdata = collect_cell_vector(V, dc)
  assemble_vector_add!(r, assembler, vecdata)

  r
end

function Gridap.Algebra.residual!(
  r::AbstractVector, odeop::DAEODEOpFromTFEOp{<:AbstractLinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache; add::Bool=false
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(odeop.tfeop)

  !add && fill!(r, zero(eltype(r)))

  yh = odeopcache.diagnostics
  yh0 = odeopcache.diagnostics0

  # Residual
  res = Gridap.ODEs.get_res(odeop.tfeop)
  dc = res(t, (uh,yh), v, yh0)

  # Forms
  order = get_order(odeop)
  forms = Gridap.ODEs.get_forms(odeop.tfeop)
  ∂tkuh = uh
  for k in 0:order
    form = forms[k+1]
    dc = dc + form(t, ∂tkuh, v)
    if k < order
      ∂tkuh = ∂t(∂tkuh)
    end
  end

  vecdata = collect_cell_vector(V, dc)
  assemble_vector_add!(r, assembler, vecdata)

  r
end

function Gridap.Algebra.allocate_jacobian(
  odeop::DAEODEOpFromTFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  Ut = evaluate(get_trial(odeop.tfeop), nothing)
  du = get_trial_fe_basis(Ut)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(odeop.tfeop)

  yh = odeopcache.diagnostics
  yh0 = odeopcache.diagnostics0

  jacs = Gridap.ODEs.get_jacs(odeop.tfeop)
  dc = DomainContribution()
  for k in 0:get_order(odeop.tfeop)
    jac = jacs[k+1]
    dc = dc + jac(t, (uh,yh), du, v, yh0)
  end
  matdata = collect_cell_matrix(Ut, V, dc)
  allocate_matrix(assembler, matdata)
end

function Gridap.ODEs.jacobian_add!(
  J::AbstractMatrix, odeop::DAEODEOpFromTFEOp,
  t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  Ut = evaluate(get_trial(odeop.tfeop), nothing)
  du = get_trial_fe_basis(Ut)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(odeop.tfeop)

  yh = odeopcache.diagnostics
  yh0 = odeopcache.diagnostics0

  jacs = Gridap.ODEs.get_jacs(odeop.tfeop)
  dc = DomainContribution()
  for k in 0:get_order(odeop)
    w = ws[k+1]
    iszero(w) && continue
    jac = jacs[k+1]
    dc = dc + w * jac(t, (uh,yh), du, v, yh0)
  end

  if num_domains(dc) > 0
    matdata = collect_cell_matrix(Ut, V, dc)
    assemble_matrix_add!(J, assembler, matdata)
  end

  J
end

function Gridap.ODEs.jacobian_add!(
  J::AbstractMatrix, odeop::DAEODEOpFromTFEOp{<:AbstractQuasilinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  Ut = evaluate(get_trial(odeop.tfeop), nothing)
  du = get_trial_fe_basis(Ut)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(odeop.tfeop)

  yh = odeopcache.diagnostics
  yh0 = odeopcache.diagnostics0

  order = get_order(odeop)
  jacs = Gridap.ODEs.get_jacs(odeop.tfeop)
  dc = DomainContribution()
  for k in 0:order-1
    w = ws[k+1]
    iszero(w) && continue
    jac = jacs[k+1]
    dc = dc + w * jac(t, (uh,yh), du, v, yh0)
  end

  # Special case for the mass matrix
  k = order
  w = ws[k+1]
  if !iszero(w)
    if Gridap.ODEs.is_form_constant(odeop, k)
      axpy_entries!(w, odeopcache.const_forms[1], J)
    else
      jac = jacs[k+1]
      dc = dc + w * jac(t, (uh,yh), du, v, yh0)
    end
  end

  if num_domains(dc) > 0
    matdata = collect_cell_matrix(Ut, V, dc)
    assemble_matrix_add!(J, assembler, matdata)
  end

  J
end

function Gridap.ODEs.jacobian_add!(
  J::AbstractMatrix, odeop::DAEODEOpFromTFEOp{<:AbstractLinearODE},
  t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
  odeopcache
)
  uh = _make_uh_from_us(odeop, us, odeopcache.Us)
  Ut = evaluate(get_trial(odeop.tfeop), nothing)
  du = get_trial_fe_basis(Ut)
  V = get_test(odeop.tfeop)
  v = get_fe_basis(V)
  assembler = Gridap.ODEs.get_assembler(odeop.tfeop)

  yh = odeopcache.diagnostics
  yh0 = odeopcache.diagnostics0

  jacs = Gridap.ODEs.get_jacs(odeop.tfeop)
  dc = DomainContribution()
  for k in 0:get_order(odeop)
    w = ws[k+1]
    iszero(w) && continue
    if Gridap.ODEs.is_form_constant(odeop, k)
      axpy_entries!(w, odeopcache.const_forms[k+1], J)
    else
      jac = jacs[k+1]
      dc = dc + w * jac(t, (uh,yh), du, v, yh0)
    end
  end

  if num_domains(dc) > 0
    matdata = collect_cell_matrix(Ut, V, dc)
    assemble_matrix_add!(J, assembler, matdata)
  end

  J
end
