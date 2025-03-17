using Symbolics

@variables α β r

Dα = Differential(α)
Dβ = Differential(β)

function mu(α,β)
  tan(β) / ( ( 1 + (tan(α))^2 + tan(β)^2  )^(1/2)     )
end

M = simplify(mu(α,β);         expand=true, threaded=false, thread_subtree_cutoff=100, rewriter=nothing )
M2 = simplify( (mu(α,β))^2;   expand=true, threaded=false, thread_subtree_cutoff=100, rewriter=nothing  )

T1 = simplify( (1-M2)^(0.5);   expand=true, threaded=false, thread_subtree_cutoff=100, rewriter=nothing )
T2 = simplify( (1-M2)^(-0.5);  expand=true, threaded=false, thread_subtree_cutoff=100, rewriter=nothing )

dMdα = simplify( expand_derivatives( Dα(M) );   expand=true, threaded=false, thread_subtree_cutoff=100, rewriter=nothing )
dMdβ = expand_derivatives( Dβ(M) )

dXdα =  -r*sin(α)*T1 - r*cos(α)*T2*M*dMdα
dXdβ = - r*cos(α)*T2*M*dMdβ

dYdα =  r*cos(α)*T1 - r*sin(α)*T2*M*dMdα
dYdβ = - r*sin(α)*T2*M*dMdβ

dZdα = r*dMdα
dZdβ = r*dMdβ


ddα = [dXdα, dYdα, dZdα]
ddβ = [dXdβ, dYdβ, dZdβ]

g11 =  ddα'*ddα
g12 =  ddα'*ddβ
g21 =  ddβ'*ddα
g22 =  ddβ'*ddβ

g = [g11 g12
    g21 g22]

ex1, ex2 = build_function(g,[α,β] )
write("metric.jl", string(ex2))
