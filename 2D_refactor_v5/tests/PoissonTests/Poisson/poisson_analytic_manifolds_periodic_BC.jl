"""
Manufacture solutions for Possion problem on doubly periodic parametric spaces
associated to a manifold (i.e. cirlce)
  - This work is only applicable in 1D to the circle with
    ϕ(x)=(r sin(x), r cos(x)) for x ∈ Ω = [0,2π] doubly periodic.
  - In 2D for a sphere, the standard chart is not an atalas, and the inverse of
    the metric is not well defined.

Solve Δu = -f on Ω, where
  - Δ is the Laplace-Beltrami operator,
  - Ω is the doubly periodic parametric space of manifold M
  - f = -surf_lap(u_ex) is a manufactured rhs

The conditions on u_ex are:
1. Domain condition: u_ex is periodic functions on Ω
2. Compatibility condition: u_ex is divergence free,
    - ∫ surf_div(surf_grad(u_ex)) = ∫surf_lap(u_ex) = 0
3. Zero mean condition:  u_ex satisfies ∫u_ex = 0

Based on previous work in flat spaces, use solve the mixed problems and enforce
∫u = ∫uex = ∫Δuex = 0 via lagrange multipler
"""


using Gridap, Gridap.Geometry
using Plots, LaTeXStrings
using LinearAlgebra
using DrWatson
using Test

include("../../../src/initialise.jl")
include("../poisson_helpers.jl")



p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns


### 1D circle of radius r
r = 2

metric_func(x) = TensorValue{1}(r^2)

u1(x) = x[1]*(2*π - x[1])
u2(x) = cos(x[1]) + 1
function u3(x)
  if x[1] < π
    return x[1]*(π-x[1])
  else
    return (x[1]-π)*(x[1]-2*π)
  end
end
function u4(x)
  if x[1] < π
    return x[1]*(π-x[1]) + x[1]*(2*π - x[1])
  else
    return (x[1]-π)*(x[1]-2*π) + x[1]*(2*π - x[1])
  end
end
u5(x) = sin(x[1]) + x[1]*(2π-x[1])

uex_funcs = Dict{Symbol,Any}()
uex_funcs[:u1] = u1
uex_funcs[:u2] = u2
uex_funcs[:u3] = u3
uex_funcs[:u4] = u4
uex_funcs[:u5] = u5



errs = []
errs_g = []
for (key, val) in uex_funcs
  for n in ns
    e, eg = solve_poisson_manifold_periodic((0,2π),(n,),p,degree,val,metric_func)
    push!(errs,e)
    push!(errs_g,eg)
  end
end

leginf = map(x->string(x),collect(keys(uex_funcs)))

# 1D circle
plot()
plot_error(ns,errs;leginf=leginf,ls=fill(:solid,length(uex_funcs)))
plot_error(ns,errs_g;ls=fill(:dash,length(uex_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)",legend=:bottomright)
plot!(ns,1e2dx.^4,lw=2,c=:blue,label="dx^4")
plot!(ns,1e1dx.^6,lw=2,c=:black,label="dx^6")
savefig(plotsdir()*"/poisson_convergence_circle_r$(Int(r))")
