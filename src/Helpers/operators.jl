

J(m::ForwardMap2Dor3D) = x -> J(m,x)
Jt(m::ForwardMap2Dor3D) = x -> Jt(m,x)
metric(m::ForwardMap2Dor3D) = x -> metric(m,x)
inv_metric(m::ForwardMap2Dor3D) = x -> inv_metric(m,x)
detg(m::ForwardMap2Dor3D)  = x -> detg(m,x)
sqrtg(m::ForwardMap2Dor3D)  = x -> sqrtg(m,x)
grad_meas(m::ForwardMap2Dor3D) = x -> gradient(sqrtg(m))(x)
covariant_basis(m::ForwardMap2Dor3D) = x -> J(m,x)
forward_jacobian(m::ForwardMap2Dor3D) = x -> J(m,x)
forward_pinv_jacobian(m::ForwardMap2Dor3D) = x -> pinvJ(J(m,x))

function pinvJ(J::MultiValue{Tuple{D1,D2}}) where {D1,D2}
  @check D2 < D1 ## J = 3x2
  Jt = transpose(J)
  inv(Jt⋅J)⋅Jt
end

function pinvJ(J::MultiValue{Tuple{D,D}}) where D
  inv(J)
end

####### surface laplacin
surflap(f::Function) = m -> surflap(f,m)
surflap(f::Function,m::ForwardMap2Dor3D) = αβ -> 1/sqrtg(m,αβ) * ( divergence(W(f,m))(αβ) )
W(f::Function,m::ForwardMap2Dor3D) = αβ ->  sqrtg(m,αβ)*( inv_metric(m,αβ) ⋅ gradient(f(m))(αβ) )

####### contra of sgrad
contr_gradf(f::Function) = m -> contr_gradf(f,m)
contr_gradf(f::Function,m::ForwardMap2Dor3D) = αβ -> inv_metric(m,αβ) ⋅ gradient(f(m))(αβ)

####### sgrad
sgrad(f::Function) = m -> sgrad(f,m)
sgrad(f::Function,m::ForwardMap2Dor3D) = αβ -> J(m,αβ) ⋅ (inv_metric(m,αβ) ⋅ gradient(f(m))(αβ) )

####### surf div
_svid(vec::Function,m::ForwardMap2Dor3D) = αβ ->  sqrtg(m,αβ)*( vec(m)(αβ))
surfdiv(vec::Function) = m -> surfdiv(vec,m)
surfdiv(f::Function,m::ForwardMap2Dor3D) = αβ -> 1/sqrtg(m,αβ) * ( divergence(_sdiv(f,m))(αβ) )

##### g⋅n for 3D thick sphere
g_star(f::Function) = m -> g_star(f,m)
g_star(f::Function,m::ForwardMap2Dor3D) = αβ -> metric(m,αβ) ⋅ VectorValue(1.0,0.0,0.0)
