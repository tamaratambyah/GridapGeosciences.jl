metric(p::Int) = αβ -> metric(p,αβ)
inv_metric(p::Int) = αβ ->  inv_metric(p,αβ)
sqrtg(p::Int) = αβ -> sqrtg(p,αβ)
grad_meas(p::Int) = αβ -> gradient(sqrtg(p))(αβ)

####### surface laplacin
W(f::Function,p::Int) = αβ ->  sqrtg(p,αβ)*( inv_metric(p,αβ) ⋅ gradient(f(p))(αβ))
surflap(f::Function,p::Int) = αβ -> 1/sqrtg(p,αβ) * ( divergence(W(f,p))(αβ) )
surflap(f::Function) = p -> surflap(f,p)

####### contra of sgrad
contr_gradf(f::Function,p::Int) = αβ -> inv_metric(p,αβ) ⋅ gradient(f(p))(αβ)
contr_gradf(f::Function) = p -> contr_gradf(f,p)

####### sgrad
sgrad(f::Function,p::Int) =  αβ -> J(p,αβ) ⋅ (inv_metric(p,αβ) ⋅ gradient(f(p))(αβ))
sgrad(f::Function) = p -> sgrad(f,p)

####### surf div
_sdiv(vec::Function,p) = αβ ->  sqrtg(p,αβ)*( vec(p)(αβ))
surfdiv(vec::Function,p::Int) = αβ -> 1/sqrtg(p,αβ) * ( divergence(_sdiv(vec,p))(αβ) )
surfdiv(vec::Function) = p -> surfdiv(vec,p)
