"""
Manufacture solutions for Possion problem on a flat domain with periodic BC.
This test is a precursor to manufacturing solutions on the sphere.

Solve Δu = -f on Ω = [0,1]^2 doubly periodic.
The weak form is:
  ∫ (∇u)⋅(∇v) dΩ = ∫ fv dΩ
where f = -Δ(u_ex) is the rhs data, and u_ex is the exact solution.

The conditions on u_ex are:
1. Domain condition: u_ex is periodic functions on [0,1]^2.
2. Compatibility condition: u_ex is divergence free, ∫ ∇⋅(∇u_ex) = ∫Δu_ex = 0
3. Zero mean condition:  u_ex satisfies ∫u_ex = 0

Examples are:
  u_ex = x(1-x) for p = 2 polynomial
      * this function is zero on x=0, x=1, and non-zero (periodic) on y=0, y=1
      * this function is NOT zero-mean, and is NOT  compatible
  u_ex = cos(2πx)
      * this function is periodic and zero-mean, and ∫Δu_ex = 0
      * this function is NOT in the FE space
  u_ex = ( x(0.5-x)     ; 0.0 < x < 0.5
           (x-0.5)(x-1) ; 0.5 < x < 1.0    )
      * this function is periodic and zero-mean, and ∫Δu_ex = 0
      * this function is in the FE space --> manufactured solutions

To stop the error blowing up, need the enforce zero mean constraint in the FE
space. This removes the zero eigenvalue which was making the linear system
singular.

Consider the following methods:
  Method 1:
    - enforce ∫u = 0 in the FE space via "FESpace(... ;constraint=zeromean)"
  Method 2:
    - enforce ∫u = 0 via lagrange multipler
    - use kwarg ;lagrange=true
  Method 3:
    - enforce ∫u = ∫uex = 0 via lagrange multipler
    - use kwarg ;lagrange=true, uzeromean=true
  Method 4: Mixed problem
    - enforce ∫u = ∫uex = ∫Δuex = 0 via lagrange multipler

Find convergence only for compatible functions (i.e. trig functions)
Find convergence only for zeromean functions (i.e. trig functions)

To manufacture solutions in FE space, consider a piecewise polynomial that has
zero mean (i.e. 'looks like a trig function')
"""

using Gridap, Gridap.Geometry
using Plots, LaTeXStrings
using LinearAlgebra
using DrWatson
dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

include("../poisson_helpers.jl")


p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns

################################################################################
#### Analytic functions
################################################################################
u1(x) = 6*x[1]*(1-x[1])
u2(x) = cos(2*π*x[1]) + 1
function u3(x)
  if x[1] < 0.5
    return x[1]*(0.5-x[1]) + 1.0
  else
    return (x[1]-0.5)*(x[1]-1.0) + 1.0
  end
end
function u4(x)
  if x[1] < 0.5
    return x[1]*(0.5-x[1]) + x[1]*(1-x[1])
  else
    return (x[1]-0.5)*(x[1]-1.0) + x[1]*(1-x[1])
  end
end
u5(x) = sin(2*π*x[1]) + x[1]*(1-x[1])

uex_funcs = Dict{Symbol,Any}()
uex_funcs[:u1] = u1
uex_funcs[:u2] = u2
uex_funcs[:u3] = u3
uex_funcs[:u4] = u4
uex_funcs[:u5] = u5


################################################################################
#### Convergence in 1D
################################################################################
method1 = []
method2 = []
method3 = []
method4 = []
for (key, val) in uex_funcs
  for n in collect(ns)
    println(n)

    e1 = solve_poisson_periodic((0,1),(n,),p,degree,val)
    e2 = solve_poisson_periodic((0,1),(n,),p,degree,val;lagrange=true)
    e3 = solve_poisson_periodic((0,1),(n,),p,degree,val;lagrange=true,uzeromean=true)
    # if n != 64
      e4 = solve_poisson_periodic_dual_form((0,1,),(n,),p,degree,val)
      push!(method4,e4)
    # end
    push!(method1,e1)
    push!(method2,e2)
    push!(method3,e3)
  end
end


leginf = map(x->string(x),collect(keys(uex_funcs)))

# method 1:
plot()
plot_error(ns,method1;leginf=leginf,ls=fill(:solid,length(uex_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)")
savefig(plotsdir()*"/poisson_convergence_periodic_method1_1D")

# method 2:
plot()
plot_error(ns,method2;leginf=leginf,ls=fill(:dashdotdot,length(uex_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)")
savefig(plotsdir()*"/poisson_convergence_periodic_method2_1D")

# method 3:
plot()
plot_error(ns,method3;leginf=leginf,ls=fill(:dot,length(uex_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)")
plot!(ns,1e1dx.^6,lw=2,c=:black,label="dx^6")
savefig(plotsdir()*"/poisson_convergence_periodic_method3_1D")

# method 4:
plot()
plot_error(ns,method4;leginf=leginf,ls=fill(:dash,length(uex_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)")
plot!(ns,1e1dx.^6,lw=2,c=:black,label="dx^6")
plot!(ns,5e-6dx.^4,lw=2,c=:blue,label="dx^4")
savefig(plotsdir()*"/poisson_convergence_periodic_method4_1D")




################################################################################
#### Convergence in 2D
################################################################################
method1 = []
method2 = []
method3 = []
method4 = []
for (key, val) in uex_funcs
  for n in collect(ns)
    println(n)

    e1 = solve_poisson_periodic((0,1,0,1),(n,n),p,degree,val)
    e2 = solve_poisson_periodic((0,1,0,1),(n,n),p,degree,val;lagrange=true)
    e3 = solve_poisson_periodic((0,1,0,1),(n,n),p,degree,val;lagrange=true,uzeromean=true)
    # if n != 64
      e4 = solve_poisson_periodic_dual_form((0,1,0,1),(n,n),p,degree,val)
      push!(method4,e4)
    # end
    push!(method1,e1)
    push!(method2,e2)
    push!(method3,e3)
  end
end


leginf = map(x->string(x),collect(keys(uex_funcs)))

# method 1:
plot()
plot_error(ns,method1;leginf=leginf,ls=fill(:solid,length(uex_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)")
savefig(plotsdir()*"/poisson_convergence_periodic_method1_2D")

# method 2:
plot()
plot_error(ns,method2;leginf=leginf,ls=fill(:dashdotdot,length(uex_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)")
savefig(plotsdir()*"/poisson_convergence_periodic_method2_2D")

# method 3:
plot()
plot_error(ns,method3;leginf=leginf,ls=fill(:dot,length(uex_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)")
plot!(ns,1e1dx.^6,lw=2,c=:black,label="dx^6")
savefig(plotsdir()*"/poisson_convergence_periodic_method3_2D")

# method 4:
plot()
plot_error(ns,method4;leginf=leginf,ls=fill(:dash,length(uex_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)")
plot!(ns,1e1dx.^6,lw=2,c=:black,label="dx^6")
plot!(ns,5e-6dx.^4,lw=2,c=:blue,label="dx^4")
savefig(plotsdir()*"/poisson_convergence_periodic_method4_2D")
