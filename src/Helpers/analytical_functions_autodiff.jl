forward_maps = map(p->ForwardMap(p),collect(1:6))
jacobians_transpose = map(m->gradient(m),forward_maps)
jacobians = map(m->Operation(transpose)(gradient(m)),forward_maps)

#### Here the point x can be 2D (α,β) or 3D (γ,α,β)
#### The correct map/jacobian will be dispatched by ForwardMap(p)
J(p::Int,x) = jacobians[p](x)
Jt(p::Int,x) =  jacobians_transpose[p](x)
metric(p::Int,x) = Jt(p,x)⋅J(p,x)
inv_metric(p::Int,x) = inv(metric(p,x))
detg(p::Int,x) = det(metric(p,x))
sqrtg(p::Int,x) = sqrt(detg(p,x))




# J(p::Int,αβ) = forward_jacobian(p)(αβ)
# Jt(p::Int,αβ) =  transpose(J(p,αβ))
# metric(p::Int,αβ) = Jt(p,αβ)⋅J(p,αβ)
# inv_metric(p::Int,αβ) = inv(metric(p,αβ))
# detg(p::Int,αβ) = det(metric(p,αβ))
# sqrtg(p::Int,αβ) = sqrt(detg(p,αβ))


########## perp operator
function perp_matrix(p::Int,αβ)
  @check length(αβ) == 2 "\n perp only defined for 2D"
  m = metric(p,αβ)
  TensorValue{2,2}( -m[1,2], m[1,1], -m[2,2], m[1,2] )
end
