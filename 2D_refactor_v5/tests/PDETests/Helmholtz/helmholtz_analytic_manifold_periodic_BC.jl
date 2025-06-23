using Gridap, Gridap.Geometry
using Plots,LaTeXStrings
using LinearAlgebra
using DrWatson
using Test

include("../analytic_metrics.jl")
include("../../../src/initialise.jl")
include("../pde_helpers.jl")

function helmholtz_manifold_periodic(manifold,n,p,degree,uex)
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
  compat = sum( ∫(ucf)dΩ + ∫( surface_laplacian(ucf,m))dΩ  )
  println("Compatibility: ", compat)

  rhs = ucf + 1.0*(surface_laplacian(ucf,m))

  #### FE Problem -- no lagrange multiplers
  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
  U = TrialFESpace(V)

  poisson_biform(u,v) = ∫(u*v)dΩg -  ∫( surface_gradient(u,m)⋅gradient(v)  )dΩg
  poisson_liform(v) = ∫(  rhs*v )dΩg
  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

  uh = solve(LUSolver(),op)

  #### Compute errors
  e = l2(uh-ucf,dΩ)
  eg = l2(uh-ucf,dΩg)
  return e,eg

end



p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns

r = 2

## on a manifold: circle
manifold =:circle

u1(x) = cos(x[1])
function u2(x)
  if x[1] < π
    return x[1]*(π-x[1])
  else
    return (x[1]-π)*(x[1]-2*π)
  end
end


uex_funcs = Dict{Symbol,Any}()
uex_funcs[:u1] = u1
uex_funcs[:u2] = u2


errs = []
errs_g = []
for (key, val) in uex_funcs
  for n in ns
    e, eg = helmholtz_manifold_periodic(manifold,n,p,degree,val)
    push!(errs,e)
    push!(errs_g,eg)
  end
end

leginf = map(x->string(x),collect(keys(uex_funcs)))

plot()
plot_error(ns,errs;leginf=leginf,ls=fill(:solid,length(uex_funcs)))
plot_error(ns,errs_g;ls=fill(:dash,length(uex_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)",legend=:bottomright)
plot!(ns,1e-1dx.^6,lw=2,c=:black,label="dx^6")
savefig(plotsdir()*"/circle_convergence")



### sphere
manifold =:sphere
u1(x) = cos(x[2])
function u3(x)
  if x[1] < 0
    return -x[1]*(π/2+x[1])
  else
    return x[1]*(x[1] - π/2)
  end
end
uex_funcs = Dict{Symbol,Any}()
uex_funcs[:u1] = u1
uex_funcs[:u2] = u3

errs = []
errs_g = []
for (key, val) in uex_funcs
  for n in ns
    e, eg = helmholtz_manifold_periodic(manifold,n,p,degree,val)
    push!(errs,e)
    push!(errs_g,eg)
  end
end

leginf = map(x->string(x),collect(keys(uex_funcs)))

plot()
plot_error(ns,errs;leginf=leginf,ls=fill(:solid,length(uex_funcs)))
plot_error(ns,errs_g;ls=fill(:dash,length(uex_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)",legend=:bottomright)
plot!(ns,1e-1dx.^6,lw=2,c=:black,label="dx^6")
savefig(plotsdir()*"/sphere_convergence")

################################################################################
################################################################################
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

poisson_biform(u,v) = ∫(u*v)dΩg -  ∫( surface_gradient(u,m)⋅gradient(v)  )dΩg
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
