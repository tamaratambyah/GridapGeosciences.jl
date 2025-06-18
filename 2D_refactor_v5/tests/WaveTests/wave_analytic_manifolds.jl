"""
Manufacture solutions for Wave equation on doubly periodic, parametric spaces
associated to a manifold.

Solve u + ∇p = fᵤ
      p + ∇⋅u = fₚ
where
  - ∇,∇⋅ is the surface gradient and surface divergence operators
  - Ω is the doubly periodic, parametric space of manifold M
  - fᵤ,fₚ are manufactured rhs

The weak form in terms of surface operators is:
  ∫ u⋅v dΩg - ∫ p* wave_div() dΩ = ∫ vfᵤ dΩg
  ∫ pq dΩg + ∫ q* surf_div(u) dΩg = ∫ vfₚ dΩg
where
  - ∫ is in the parametric space
  - surf_grad, surf_div, dΩg account for the metric
  - wave_div is a new operator
  - dΩ is the standard measure of the parametric space

Consider various mappings:
  1D interval -> 2D circle
    * φ(x)  = (r sin(x), r cos(x))
    * x ∈ Ω = [0, 2π] periodic
    * u_ex  = x(2π - x)
    * p_ex  = 1 + 0.01x(2π-x)

  2D rectangle -> 3D sphere
    * φ(u,v)  = (r cos(u) cos(v), r cos(u) sin(v), r sin(u)  )
    * u,v ∈ Ω = [-π/2, π/2] × [0, 2π] doubly periodic
    * u_ex    = (0, v(2π - v))
    * p_ex    = 1 + 0.01(π/2 + u)(π/2 - u)
"""

using Gridap, Gridap.Geometry
using Plots, LaTeXStrings
using LinearAlgebra
using DrWatson

include("../../src/initialise.jl")
include("wave_helpers.jl")
include("../PoissonTests/poisson_helpers.jl")



p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns

radii = [0.1, 0.5, 1.0, 2.0, 4.0]

################################################################################
#### 1D tests: circle
################################################################################
# uex(x) = VectorValue(x[1]*(2π-x[1]))
# pex(x) = 1.0 + 0.01*x[1]*(2π-x[1])

uex(x) = VectorValue(cos(x[1]))
pex(x) = 1.0 + 0.1*sin(x[1])

errs_u = []
errs_ug = []
errs_p = []
errs_pg = []

for r in collect(radii)
  metric_func(x) = TensorValue{1}(r^2)
  for n in collect(ns)
    println(n)

    eu,ep,eu_g,ep_g = solve_wave_manifold((0,2π),(n, ),p,degree,metric_func,uex,pex)

    push!(errs_u,eu)
    push!(errs_ug,eu_g)
    push!(errs_p,ep)
    push!(errs_pg,ep_g)
  end
end

leginf = map(x->"r = $x",collect(radii))

plot()
plot_error(ns,errs_u;leginf=leginf,ls=fill(:solid,length(leginf)))
plot_error(ns,errs_ug;ls=fill(:dash,length(leginf)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)")
plot!(ns,1e1dx.^8,lw=2,c=:black,label="dx^8")
savefig(plotsdir()*"/wave_convergence_u_1D_trig")

plot()
plot_error(ns,errs_p;leginf=leginf,ls=fill(:solid,length(leginf)))
plot_error(ns,errs_pg;ls=fill(:dash,length(leginf)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(p - ph)")
plot!(ns,dx.^6,lw=2,c=:black,label="dx^6")
savefig(plotsdir()*"/wave_convergence_p_1D_trig")

################################################################################
#### 2D tests: sphere
################################################################################
# uex(x) = VectorValue(0.0,x[2]*(2π-x[2]))
# pex(x) = 1.0 +  0.01*(π/2+x[1])*(π/2-x[1])

uex(x) = VectorValue(cos(x[1]),cos(x[2]))
pex(x) = 1.0 +  0.1*cos(x[1])

errs_u = []
errs_ug = []
errs_p = []
errs_pg = []
for r in collect(radii)
  metric_func(x) = TensorValue{2,2}(r^2 ,0.0, 0.0, r^2*(cos(x[1]))^2 )
  for n in collect(ns)
    println(n)

    eu,ep,eu_g,ep_g = solve_wave_manifold((-π/2,π/2, 0,2π),(n,n),p,degree,metric_func,uex,pex)

    push!(errs_u,eu)
    push!(errs_ug,eu_g)
    push!(errs_p,ep)
    push!(errs_pg,ep_g)
  end
end

leginf = map(x->"r = $x",collect(radii))

plot()
plot_error(ns,errs_u;leginf=leginf,ls=fill(:solid,length(leginf)))
plot_error(ns,errs_ug;ls=fill(:dash,length(leginf)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)")
plot!(ns,1e2dx.^8,lw=2,c=:black,label="dx^8")
savefig(plotsdir()*"/wave_convergence_u_2D_trig")

plot()
plot_error(ns,errs_p;leginf=leginf,ls=fill(:solid,length(leginf)))
plot_error(ns,errs_pg;ls=fill(:dash,length(leginf)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(p - ph)")
plot!(ns,dx.^6,lw=2,c=:black,label="dx^6")
savefig(plotsdir()*"/wave_convergence_p_2D_trig")
