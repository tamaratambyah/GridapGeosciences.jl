"""
Manufacture solutions for the helmholtz problem on doubly periodic parametric spaces
associated to a manifold (i.e. cirlce/sphere)

Solve u + Δu = f on Ω, where
  - Δ is the Laplace-Beltrami operator,
  - Ω is the doubly periodic parametric space of manifold M
  - f = uex + surf_lap(u_ex) is a manufactured rhs, such that ∫f = 0

The conditions on u_ex are:
1. Domain condition: u_ex is periodic functions on Ω
2. Compatibility condition: u_ex is divergence free,
    - ∫ surf_div(surf_grad(u_ex)) = ∫surf_lap(u_ex) = 0
3. Zero mean condition:  u_ex satisfies ∫u_ex = 0

The condition on f is ∫f = ∫(uex+Δuex) = 0, to ensure manufactured solutions.
This problem has no kernel (i.e. no zero eigvals)
"""

using Gridap, Gridap.Geometry
using Plots,LaTeXStrings
using LinearAlgebra
using DrWatson
using Test


include("../analytic_metrics.jl")
include("../../../src/initialise.jl")
include("../pde_helpers.jl")
include("HelmholtzSolvers.jl")


p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]


r = 2

## on a manifold: circle
errs = []
errs_g = []

errs_K = []
errs_gK = []
for (key, val) in uex_periodic_funcs
  for n in ns
    e, eg = helmholtz_manifold_periodic(:circle,n,p,degree,val)
    eK, egK = helmholtz_manifold_periodic_symmetricK(:circle,n,p,degree,val)
    push!(errs,e)
    push!(errs_g,eg)
    push!(errs_K,eK)
    push!(errs_gK,egK)
  end
end

leginf = map(x->string(x),collect(keys(uex_periodic_funcs)))
dx = (2π ./ ns )
plot()
plot_error(ns,errs;leginf=leginf,ls=fill(:solid,length(uex_periodic_funcs)))
plot_error(ns,errs_g;ls=fill(:dash,length(uex_periodic_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)",legend=:bottomright)
plot!(ns,1e-1dx.^6,lw=2,c=:black,label="dx^6")
savefig(plotsdir()*"/circle_convergence")

plot()
plot_error(ns,errs_K;leginf=leginf,ls=fill(:solid,length(uex_periodic_funcs)))
plot_error(ns,errs_gK;ls=fill(:dash,length(uex_periodic_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)",legend=:bottomright)
plot!(ns,1e-1dx.^6,lw=2,c=:black,label="dx^6")
savefig(plotsdir()*"/circle_convergence_symmetricK")


### sphere
errs = []
errs_g = []

errs_K = []
errs_gK = []
for (key, val) in uex_periodic_funcs
  for n in ns
    e, eg = helmholtz_manifold_periodic(:sphere,n,p,degree,val)
    eK, egK = helmholtz_manifold_periodic_symmetricK(:circle,n,p,degree,val)
    push!(errs,e)
    push!(errs_g,eg)
    push!(errs_K,eK)
    push!(errs_gK,egK)
  end
end

leginf = map(x->string(x),collect(keys(uex_periodic_funcs)))
dx = (2π ./ ns ) .* (π ./ ns)
plot()
plot_error(ns,errs;leginf=leginf,ls=fill(:solid,length(uex_periodic_funcs)))
plot_error(ns,errs_g;ls=fill(:dash,length(uex_periodic_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)",legend=:bottomright)
plot!(ns,1e-1dx.^3,lw=2,c=:black,label="dx^3")
savefig(plotsdir()*"/sphere_convergence")

plot()
plot_error(ns,errs_K;leginf=leginf,ls=fill(:solid,length(uex_periodic_funcs)))
plot_error(ns,errs_gK;ls=fill(:dash,length(uex_periodic_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)",legend=:bottomright)
plot!(ns,1e-1dx.^6,lw=2,c=:black,label="dx^6")
savefig(plotsdir()*"/sphere_convergence_symmetricK")




################################################################################
##### low level
################################################################################
n = 4
# uex(x) = cos(x[1])
function uex(x)
  if x[1] < π
    return x[1]*(π-x[1])
  else
    return (x[1]-π)*(x[1]-2*π)
  end
end
manifold=:sphere
_metric_func = metrics[manifold]
domain = domains[manifold]


d = Int(length(domain)/2)

model = CartesianDiscreteModel(domain, ntuple(x->n,d), isperiodic=ntuple(x->true,d))
Ω = Triangulation(model)
m = Metric(_metric_func,Ω)

dΩ = Measure(Ω,degree)
dΩg =  Measure(m,Ω,degree)


ucf = CellField(uex,Ω)


# check compatibility
sum( ∫(ucf)dΩ + ∫( surface_laplacian(ucf,m))dΩ  )

writevtk(Ω,dir*"/poisson",cellfields=["u"=>uex],append=false)


_rhs = ucf + 1.0*(surface_laplacian(ucf,m))

#### FE Problem -- no lagrange multiplers
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

D(x) = TensorValue{3,2}(-r*sin(x[1])*cos(x[2]),r*cos(x[2])*cos(x[1]),0,
                        -r*cos(x[1])*sin(x[2]),-r*sin(x[2])*sin(x[1]),r*cos(x[2]) )

stiffnes(u,v) =    ∫( (D ⋅ ((m.inv_metric⋅ gradient(v)) ) ⋅( D⋅(m.inv_metric⋅ gradient(u))) )   )dΩg
K = assemble_matrix(stiffnes,U,V)
issymmetric(Array(K))

stiffness2(u,v) =  ∫( surface_gradient(v,m)⋅gradient(u)  )dΩg
K2 = assemble_matrix(stiffness2,U,V)
issymmetric(Array(K2))
Array(K) == Array(K2)


poisson_biform(u,v) = ∫(u*v)dΩg - stiffnes(u,v)#   ∫( surface_gradient(u,m)⋅gradient(v)  )dΩg
poisson_liform(v) = ∫(  _rhs*v )dΩg
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

A = get_matrix(op)
b = get_vector(op)

evals = eigvals(Array(A))
sum(A*(3*ones(size(b))))
sum(b)

uh = solve(LUSolver(),op)
sum( ∫(uh )dΩ  )

#### Compute errors
e = uh-uex
sum(∫(e⊙e)dΩ)
sum(∫(e⊙e)dΩg)

writevtk(Ω,dir*"/poisson",
        cellfields=["u"=>uex,"uh"=>uh,"e"=>e],append=false)
