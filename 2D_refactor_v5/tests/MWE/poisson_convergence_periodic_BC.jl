"""
Manufacture solutions for Possion problem on a flat domain with periodic BC.
This test is a precursor to manufacturing solutions on the sphere.

Solve Δu = -f on Ω = [0,1]^2 doubly periodic.
The weak form is:
  ∫ (∇u)⋅(∇v) dΩ = ∫ fv dΩ
where f = -Δ(u_ex) is the rhs data.

Since the domain is doubly periodic, there are no boundary integrals in the weak
form. In addition, consider u_ex as periodic functions on [0,1]^2.
Examples are:
  u_ex = x(1-x) for p = 2 polynomial
      * this function is zero on x=0, x=1, and non-zero (periodic) on y=0, y=1
      * this function is NOT zero-mean
  u_ex = cos(2πx)
      * this function is periodic and zero-mean
      * this function is NOT in the FE space
  u_ex = ( x(0.5-x)     ; 0.0 < x < 0.5
           (x-0.5)(x-1) ; 0.5 < x < 1.0    )
      * this function is periodic and zero-mean
      * this function is in the FE space (use to manufacture solutions)

To stop the error blowing up, need the enforce zero mean constraint in the FE
space. This removes the zero eigenvalue which was making the linear system
singular.
Test 2 ways of enforcing zero mean
  1. in the FE space via "FESpace(... ;constraint=zeromean")
  2. algebraically via lagrange multipler
Find that both methods give the same convergence results (expected)
Find convergence only for zeromean functions (i.e. trig functions)
  - Trying to force zero mean to a non zeromean function does not yield good
    convergence results.
To manufacture solutions in FE space, consider a piecewise polynomial that has
zero mean (i.e. 'looks like a trig function')
"""

using Gridap, Gridap.Geometry
using Plots
using LinearAlgebra
using DrWatson
dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

function solve_poisson_periodic(domain,n,p,degree,uex)

  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,n,isperiodic=(true,true)))

  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  ## assess zero mean of analytic function
  println("Zero mean is: " , sum(∫(uex)dΩ))

  #### FE Problem

  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1,
                        constraint=:zeromean)
  U = TrialFESpace(V)

  rhs(x) = -1.0*(laplacian(uex)(x))

  poisson_biform(u,v) = ∫( gradient(u)⋅gradient(v) )dΩ
  poisson_liform(v) =  ∫( rhs*v )dΩ

  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
  uh = solve(LUSolver(),op)

  e = sum(∫((uh-uex)⊙(uh-uex))dΩ)

  println("Errors: ", e)
  return e
end

function solve_poisson_periodic_lagrange(domain,n,p,degree,uex)

  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,n,isperiodic=(true,true)))

  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  ## assess zero mean of analytic function
  println("Zero mean is: " , sum(∫(uex)dΩ))

  #### FE Problem

  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
  U = TrialFESpace(V)
  Λ = ConstantFESpace(model)
  M = TrialFESpace(Λ)
  X = MultiFieldFESpace([U,M])
  Y = MultiFieldFESpace([V,Λ])

  rhs(x) = -1.0*(laplacian(uex)(x))

  poisson_biform((u,μ),(v,λ)) = ∫( gradient(u)⋅gradient(v)  )dΩ  + ∫(v*μ)dΩ + ∫(λ*u)dΩ
  poisson_liform((v,λ)) = ∫( rhs*v )dΩ + ∫(λ*uex)dΩ

  op = AffineFEOperator(poisson_biform,poisson_liform,X,Y)
  uh,μh = solve(LUSolver(),op)

  e = sum(∫((uh-uex)⊙(uh-uex))dΩ)

  println("Errors: ", e)
  return e
end

function plot_error(ns,errs,leginf;
    ls=[:solid, :dash, :dot, :dashdot, :dashdotdot],
    colors = palette(:tab10),
    markers = [:circle, :rect, :diamond, :utriangle, :cross, :xcross] )

  nsims = Int(length(errs)/length(ns))
  for i in 1:nsims
    idx1 = 1 + (i-1)*length(ns)
    idx2 = (i)*length(ns)
    plot!(ns,
          errs[idx1:idx2],
          lw=3,
          c=colors[i],ls=ls[i], markershape=markers[i],
          label=leginf[i])
  end

end

#### Analytic solution with zero mean
p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns

dd = Dict(
          # "x(1-x)" => ( u0(x) = x[1]*(1-x[1]) ),
          "cos(2πx)" => ( u1(x) = cos(2*π*x[1])  ),
          "sin(2πx)" => ( u2(x) =  sin(2*π*x[1]) ),
          "cos(2πx)cos(2πy)" => ( u3(x) =  cos(2*π*x[1])*cos(2*π*x[2]) )
 )


errs = []
errs_lagrange = []
for (key, val) in dd
  for n in collect(ns)
    e = solve_poisson_periodic((0,1,0,1),(n,n),p,degree,val)
    eg = solve_poisson_periodic_lagrange((0,1,0,1),(n,n),p,degree,val)
    push!(errs,e)
    push!(errs_lagrange,eg)
  end
end

plot()
plot_error(ns,errs,collect(keys(dd));ls=fill(:solid,length(dd)))
plot_error(ns,errs_lagrange,collect(keys(dd));ls=fill(:dashdotdot,length(dd)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - u_h)"
)
plot!(ns,5e-1dx.^6,lw=2,c=:black,ls=:dash,label="dx^6")
savefig(plotsdir()*"/poisson_convergence_periodic")

################# piecewise zero mean polynomial function
function uex(x)
  if x[1] < 0.5
    return x[1]*(0.5-x[1])
  else
    return (x[1]-0.5)*(x[1]-1.0)
  end
end

errs = []
errs_lagrange = []
for n in collect(ns)
  e = solve_poisson_periodic((0,1,0,1),(n,n),p,degree,uex)
  eg = solve_poisson_periodic_lagrange((0,1,0,1),(n,n),p,degree,uex)
  push!(errs,e)
  push!(errs_lagrange,eg)
end


plot()
plot_error(ns,errs,["constraint"])
plot_error(ns,errs_lagrange,["lagrange"];ls=[:dash])
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - u_h)"
)
savefig(plotsdir()*"/poisson_manufactured_periodic")



################# force zero-mean into non zero mean analytic function
u(x) = x[1]*(1-x[1])

model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(16,16),isperiodic=(true,true)))
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)
u_mean = sum(∫( u )dΩ)
uex(x) = u(x) - u_mean

errs = []
errs_lagrange = []
for n in collect(ns)
  e = solve_poisson_periodic((0,1,0,1),(n,n),p,degree,uex)
  eg = solve_poisson_periodic_lagrange((0,1,0,1),(n,n),p,degree,uex)
  push!(errs,e)
  push!(errs_lagrange,eg)
end
plot()
plot_error(ns,errs,["constraint"])
plot_error(ns,errs_lagrange,["lagrange"];ls=[:dash])
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - u_h)"
)
savefig(plotsdir()*"/poisson_convergence_periodic_poly_force")





################# low level
# u(x) =  x[1]*(1-x[1]) #cos(2*π*x[1])*cos(2*π*x[2])
# _u(x) = u(x) - sum(∫( u )dΩ)
degree = 2
p = 2
function u(x)
  if x[1] < 0.5
    return x[1]*(0.5-x[1])
  else
    return (x[1]-0.5)*(x[1]-1.0)
  end
end


rhs(x) = -1.0*(laplacian(u)(x))

#### Domain and quadrature
degree = 2*(p+1)

model = CartesianDiscreteModel((0,1, 0,1), (12,12), isperiodic=(true,true))

Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

sum(∫(u)dΩ  )

writevtk(Ω,dir*"/poisson_manufactured_periodic_BC",
        cellfields=["u"=>u],append=false)

## assess zero mean of analytic function
sum(∫(u)dΩ)

#### FE Problem -- no lagrange multiplers
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1,constraint=:zeromean)
U = TrialFESpace(V)

poisson_biform(u,v) = ∫( gradient(u)⋅gradient(v)  )dΩ
poisson_liform(v) = ∫(  rhs*v )dΩ
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

sum( ∫(uh )dΩ  )

#### Compute errors
e = uh-u
sum(∫(e⊙e)dΩ)

writevtk(Ω,dir*"/poisson_manufactured_periodic_BC",
        cellfields=["u"=>u,"uh"=>uh,"e"=>e],append=false)
