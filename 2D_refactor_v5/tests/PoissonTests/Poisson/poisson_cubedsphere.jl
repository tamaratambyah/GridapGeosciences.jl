"""
Manufacture solutions for Possion problem on the cubed sphere mesh.

Solve Δu = -f on Ω, where
  - Δ is the Laplace-Beltrami operator,
  - Ω is the parametric space of manifold M
  - f = -surf_lap(u_ex) is a manufactured rhs

The weak form in terms of surface operators is:
  ∫ (surf_grad(u)) ⋅ (grad(v)) dΩg = ∫ v(-surf_lap(u_ex)) dΩg
where
  - ∫ is in the parametric space
  - surf_grad, surf_lap, dΩg account for the metric
  - grad is the standard flat gradient in the parametric space

The parametric space is Ω = [-π/4,π/4]^2 doubly periodic
  - need to enforce zero mean constraint
Consider piecewise u_ex that is periodic, zeromean and in FE space:
  - u_ex = ( -α*(α + π/4 ) if α < 0 ;
              α*(α - π/4 ) is α > 0 )
Find that the error is large
"""

using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields
using Plots, LaTeXStrings
include("../../../src/initialise.jl")
include("../poisson_helpers.jl")


p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns

################################################################################
#### Analytic parametric space for a single panel
################################################################################

global RADIUS = 1.0*sqrt(3.0)


function u1(αβ)
  if αβ[1] < 0.0
    return -αβ[1]*(αβ[1] + π/4)
  else
    return αβ[1]*(αβ[1] - π/4)
  end
end
u2(x) = -(x[1] + π/4  )*(x[1] - π/4)
u3(x) = cos(4*x[1])
# n=4
# e, eg = solve_poisson_manifold_periodic((-π/4,π/4, -π/4,π/4), (n,n),p,degree,u1,metric_func(cubedsphere);name=cubedsphere)


uex_funcs = Dict{Symbol,Any}()
uex_funcs[:u1] = u1
uex_funcs[:u2] = u2
uex_funcs[:u3] = u3

errs = []
errs_g = []
for (key, val) in uex_funcs
  for n in ns[1:end-1]
    e, eg = solve_poisson_manifold_periodic((-π/4,π/4, -π/4,π/4), (n,n),p,degree,val,metric_func(cubedsphere);name=cubedsphere)
    push!(errs,e)
    push!(errs_g,eg)
  end
end


leginf = map(x->string(x),collect(keys(uex_funcs)))

#
plot()
plot_error(ns[1:end-1],errs;leginf=leginf,ls=fill(:solid,length(uex_funcs)))
plot_error(ns[1:end-1],errs_g;ls=fill(:dash,length(uex_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)",legend=:bottomright)
savefig(plotsdir()*"/poisson_cubedsphere_analytic_domain_new")
