using Gridap
using GridapGeosciences
using Test

### some points for testing
function _r(x)
  a = -π/4 + (π/4 - -π/4)*rand()
  b = -π/4 + (π/4 - -π/4)*rand()
  Point(a,b)
end
pts = [_r(1) for i in 1:10]
αβ = pts[1]
p = 1


########### Forward jacobian
auto_forward_jacobian(p::Int,αβ) = transpose( gradient(forward_map_2D(p))(αβ) )
auto_forward_jacobian(p::Int) = αβ -> auto_forward_jacobian(p::Int,αβ)

auto_forward_jacobian(p,αβ) ≈ forward_jacobian(αβ,p)
auto_forward_jacobian(p)(αβ) ≈ forward_jacobian(p)(αβ)

# auto_forward_jacobian(p) = αβ -> transpose( gradient(forward_map_2D(p))(αβ) )

############### ONE VARIABLES ##################################################
# J1(αβ) =  auto_forward_jacobian(1)(αβ)
# J1t(αβ) =  transpose(J1(αβ))
# auto_metric(αβ) = J1t(αβ)⋅J1(αβ)
# auto_inv_metric(αβ) = inv(auto_metric(αβ))
# auto_detg(αβ) = det(auto_metric(αβ))
# auto_sqrtg(αβ) = sqrt(auto_detg(αβ))
# auto_grad_meas(αβ) = gradient(auto_sqrtg)(αβ)

# auto_grad_meas(αβ)

# ####### surface laplacin
# auto_W(f::Function,p::Int) = αβ ->  auto_sqrtg(αβ)*( auto_inv_metric(αβ) ⋅ gradient(f(p))(αβ))
# auto_surflap(f::Function,p::Int) = αβ -> 1/auto_sqrtg(αβ) * ( divergence(auto_W(f,p))(αβ) )
# auto_surflap(f::Function) = p -> auto_surflap(f,p)

# p=2
# auto_surflap(f_sin)(p)(αβ)
# surflap(f_sin)(p)(αβ)


############### TWO VARIABLES ##################################################

J(p::Int,αβ) = auto_forward_jacobian(p)(αβ)
Jt(p::Int,αβ) =  transpose(J(p,αβ))
auto_metric(p::Int,αβ) = Jt(p,αβ)⋅J(p,αβ)
auto_inv_metric(p::Int,αβ) = inv(auto_metric(p,αβ))
auto_detg(p::Int,αβ) = det(auto_metric(p,αβ))
auto_sqrtg(p::Int,αβ) = sqrt(auto_detg(p,αβ))
auto_sqrtg(p::Int) = αβ -> auto_sqrtg(p,αβ)
auto_grad_meas(p::Int) = αβ -> gradient(auto_sqrtg(p))(αβ)

auto_grad_meas(p)(αβ)

####### surface laplacin
auto_W(f::Function,p::Int) = αβ ->  auto_sqrtg(p,αβ)*( auto_inv_metric(p,αβ) ⋅ gradient(f(p))(αβ))
auto_surflap(f::Function,p::Int) = αβ -> 1/auto_sqrtg(p,αβ) * ( divergence(auto_W(f,p))(αβ) )
auto_surflap(f::Function) = p -> auto_surflap(f,p)

p=1
auto_surflap(f_sin)(p)(αβ)
surflap(f_sin)(p)(αβ)

####### contra of sgrad
auto_contr_gradf(f::Function,p::Int) = αβ -> auto_inv_metric(p,αβ) ⋅ gradient(f(p))(αβ)
auto_contr_gradf(f::Function) = p -> auto_contr_gradf(f,p)
auto_contr_gradf(f_sin)(p)(αβ) ≈ contr_gradf(f_sin)(p)(αβ)

####### sgrad
auto_sgrad(f::Function,p::Int) =  αβ -> J(p,αβ) ⋅ (auto_inv_metric(p,αβ) ⋅ gradient(f(p))(αβ))
auto_sgrad(f::Function) = p -> auto_sgrad(f,p)

auto_sgrad(f_sin)(p)(αβ) ≈ sgrad(f_sin)(p)(αβ)


####### surf div
_a_sdiv(vec::Function,p) = αβ ->  auto_sqrtg(p,αβ)*( vec(p)(αβ))
auto_surfdiv(vec::Function,p::Int) = αβ -> 1/auto_sqrtg(p,αβ) * ( divergence(_a_sdiv(vec,p))(αβ) )
auto_surfdiv(vec::Function) = p -> auto_surfdiv(vec,p)

auto_surfdiv(contr_gradf(f_sin))(p)(αβ) ≈ surfdiv(contr_gradf(f_sin))(p)(αβ)





import Gridap.TensorValues: MultiValue
function pinvJ(J::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  Jt = transpose(J)
  inv(Jt⋅J)⋅Jt
end
auto_forward_pinv_jacobian(p) = αβ -> pinvJ(forward_jacobian(p)(αβ))
auto_forward_pinv_jacobian(p)(αβ)  ≈ forward_pinv_jacobian(p)(αβ)

vX = panel_to_cartesian(tangent_vec(vecX))

auto_contra_v(vecX::Function,p::Int) = αβ -> auto_forward_pinv_jacobian(p)(αβ)⋅ vecX(p)(αβ)
auto_contra_v(vecX::Function) = p -> auto_contra_v(vecX,p)


contra_v(vX)(p)(αβ) ≈ auto_contra_v(vX)(p)(αβ)
