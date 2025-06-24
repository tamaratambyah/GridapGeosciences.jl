using Gridap, Gridap.Geometry
using Plots,LaTeXStrings
using LinearAlgebra
using DrWatson
using Test

include("../pde_helpers.jl")
include("HelmholtzSolvers.jl")

p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]


#### compatible analytic functions
function u1(x)
  if x[1] < 0.5
    return x[1]*(0.5-x[1])
  else
    return (x[1]-0.5)*(x[1]-1)
  end
end

u2(x) = cos(2*π*x[1])

uex_funcs = Dict{Symbol,Any}()
uex_funcs[:u1] = u1
uex_funcs[:u2] = u2

errs = []
errs_dual = []
for (key, val) in uex_funcs
  for n in ns
    e = helmholtz_periodic((0,1,0,1),(n,n),p,degree,val)
    e_d = helmholtz_dual_periodic((0,1,0,1),(n,n),p,degree,val)
    push!(errs,e)
    push!(errs_dual,e_d)
  end
end

leginf = map(x->string(x),collect(keys(uex_funcs)))
dx = 1 ./ (ns.^2)

plot()
plot_error(ns,errs;leginf=leginf,ls=fill(:solid,length(uex_funcs)))
plot_error(ns,errs_dual;ls=fill(:dash,length(uex_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)",legend=:bottomright)
plot!(ns,1e1dx.^3,lw=2,c=:black,label="dx^3")
savefig(plotsdir()*"/helmholtz_convergence")



################################################################################
################################################################################
model = CartesianDiscreteModel((0,1,0,1), (8,8), isperiodic=(true,true))
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

####### Primial form
uex = cos(2*π*x[1])

ucf = CellField(uex,Ω)
# check compatibility
( sum( ∫( ucf)dΩ) + sum(∫( laplacian(ucf))dΩ  ) )

writevtk(Ω,dir*"/poisson",cellfields=["u"=>uex],append=false)


rhs = (ucf + laplacian(ucf))

#### FE Problem -- no lagrange multiplers
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

poisson_biform(u,v) = ∫(u*v)dΩ -  ∫( gradient(u)⋅gradient(v)  )dΩ
poisson_liform(v) = ∫(  rhs*v )dΩ
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

A = get_matrix(op)
b = get_vector(op)

evals = eigvals(Array(A))
sum(A*(3*ones(size(b))))
sum(b)

uh = solve(LUSolver(),op)
# sum( ∫(uh )dΩ  )

#### Compute errors
e = uh-uex
sum(∫(e⊙e)dΩ)

writevtk(Ω,dir*"/poisson",
        cellfields=["u"=>uex,"uh"=>uh,"e"=>e],append=false)





#### Mixed form
# uex = cos(2*π*x[1])
function uex(x)
  if x[1] < 0.5
    return x[1]*(0.5-x[1])
  else
    return (x[1]-0.5)*(x[1]-1)
  end
end

ucf = CellField(uex,Ω)
# check compatibility
( sum( ∫( ucf)dΩ) + sum(∫( laplacian(ucf))dΩ  ) )

# dual form
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:L2)
U = TrialFESpace(V)

T = TestFESpace(Ω, ReferenceFE(raviart_thomas,Float64,p); conformity=:Hdiv)
S = TrialFESpace(T)

X = MultiFieldFESpace([S,U])
Y = MultiFieldFESpace([T,V])


_rhs = (ucf + laplacian(ucf))

biformX((s,u),(t,v)) = (  ∫( s⋅t)dΩ + ∫( divergence(t)*u )dΩ
                        + ∫( u*v )dΩ   + ∫( divergence(s)*v  )dΩ
                      )
liformY((t,v)) = ∫( _rhs*v )dΩ

op = AffineFEOperator(biformX,liformY,X,Y)
sh,uh = solve(LUSolver(),op)

e = uh-uex
sum(∫(e⊙e)dΩ)
