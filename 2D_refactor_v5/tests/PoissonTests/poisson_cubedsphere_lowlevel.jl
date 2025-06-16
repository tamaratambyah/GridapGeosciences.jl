using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields
using Plots, LaTeXStrings
include("../../src/initialise.jl")
include("poisson_helpers.jl")



#### Analytic solution with zero mean
p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns

function uzeroemean(αβ)
  cos(4*αβ[1])
  # if αβ[1] < 0.0
  #   return -αβ[1]*(αβ[1] + π/4)
  # else
  #   return αβ[1]*(αβ[1] - π/4)
  # end
end

### no metric
errs = []
errs_lagrange = []
for n in collect(ns)
  e = solve_poisson_periodic((-π/4,π/4, -π/4,π/4),(n,n),p,degree,uzeroemean)
  eg = solve_poisson_periodic_lagrange((-π/4,π/4, -π/4,π/4),(n,n),p,degree,uzeroemean)
  push!(errs,e)
  push!(errs_lagrange,eg)
end

plot()
plot_error(ns,errs,["zeromean"];ls=fill(:solid,1))
plot_error(ns,errs_lagrange,["lagrange"];ls=fill(:dashdotdot,1))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - u_h)"
)
plot!(ns,5e-1dx.^6,lw=2,c=:black,ls=:dash,label="dx^6")
savefig(plotsdir()*"/poisson_convergence_periodic")

### with metric -- does not converge
# global A_bump, B_bump, b_bump = bump_matrics(π/4)
global RADIUS = 1.0*sqrt(3.0)
errs = []
errs_lagrange = []
for n in collect(ns)
  e = solve_poisson_manifold((-π/4,π/4, -π/4,π/4),(n,n),p,degree,uzeroemean,metric_func(cubedsphere);isperiodic=(true,true))
  eg = solve_poisson_manifold((-π/4,π/4, -π/4,π/4),(n,n),p,degree,uzeroemean,metric_func(cubedsphere);isperiodic=(true,true))
  push!(errs,e)
  push!(errs_lagrange,eg)
end
