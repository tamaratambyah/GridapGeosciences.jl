import Gridap.Helpers: @abstractmethod
import Gridap.ODEs: LinearODE, QuasilinearODE, SemilinearODE
import Gridap.ODEs: ODEOperatorType
import Gridap.ODEs: ODEOperator

###############
# DAEODEOperator #
###############

abstract type DAEODEOperator{T} <: ODEOperator{T} end
ODEOperatorType(::DAEODEOperator{T}) where {T} = T
ODEOperatorType(::Type{<:DAEODEOperator{T}}) where {T} = T

function Gridap.Algebra.residual(
  odeop::DAEODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}},
  odeopcache
)
  r = allocate_residual(odeop, t, us, odeopcache)
  Gridap.ODEs.residual!(r, odeop, t, us, odeopcache)
  r
end




function Gridap.Algebra.jacobian!(
  J::AbstractMatrix, odeop::DAEODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
  odeopcache
)
  LinearAlgebra.fillstored!(J, zero(eltype(J)))
  Gridap.ODEs.jacobian_add!(J, odeop, t, us, ws, odeopcache)
  J
end


function Gridap.Algebra.jacobian(
  odeop::DAEODEOperator,
  t::Real, us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}},
  odeopcache
)
  J = allocate_jacobian(odeop, t, us, odeopcache)
  Gridap.ODEs.jacobian!(J, odeop, t, us, ws, odeopcache)
  J
end
