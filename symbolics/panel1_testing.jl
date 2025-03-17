using DrWatson
using Symbolics
using Latexify

# u = Symbolics.variables(:u, 1:2)
@variables α,β

@variables x

D = Differential(x)
f(x) = (tan(x))^2
expand_derivatives(D(f(x)))

g(x) = asin(x)
expand_derivatives(D(g(x)))

function gamma(α,β)
  # α,β = αβ
  [α,
   asin( ( tan(β) )/( ( 1 + (tan(α))^2 + (tan(β))^2 )^(0.5) ) ) ]
end


function sigma(θϕ)
  θ,ϕ = θϕ

  [cos(θ)*cos(ϕ),
   sin(θ)*cos(ϕ),
   sin(ϕ)]

end

@variables θ ϕ α β
B = [cos(θ)*cos(ϕ),
sin(θ)*cos(ϕ),
sin(ϕ)]

C = simplify.(substitute.(B, (Dict(θ => α, ϕ => asin( ( tan(β) )/( ( 1 + (tan(α))^2 + (tan(β))^2 )^(0.5) ) )
  ), )))

D1 = Differential.([α,β]);
jac1 = expand_derivatives.(collect(Symbolics.@arrayop (i,j) D1[j](C[i])))


jac2 = Symbolics.jacobian(C, [α,β];simplify=true)




kappa(α,β) = simplify.( sigma(gamma(α,β));
  expand=true,
threaded=false,
thread_subtree_cutoff=100,
rewriter=nothing   )

kappa(α,β)
