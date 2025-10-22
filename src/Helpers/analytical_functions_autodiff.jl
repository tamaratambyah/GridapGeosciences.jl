

J(p::Int,αβ) = forward_jacobian(p)(αβ)
Jt(p::Int,αβ) =  transpose(J(p,αβ))
metric(p::Int,αβ) = Jt(p,αβ)⋅J(p,αβ)
inv_metric(p::Int,αβ) = inv(metric(p,αβ))
detg(p::Int,αβ) = det(metric(p,αβ))
sqrtg(p::Int,αβ) = sqrt(detg(p,αβ))


########## perp operator
function perp_matrix(p::Int,αβ)
  m = metric(p,αβ)
  TensorValue{2,2}( -m[1,2], m[1,1], -m[2,2], m[1,2] )
end
