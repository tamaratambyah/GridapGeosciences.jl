using Gridap, Gridap.Geometry
using Plots
using LinearAlgebra
using DrWatson
using Test
dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


function poisson_dual_form(domain,partition,p,degree,uex)

  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,partition,isperiodic=ntuple(x->true,length(partition))))

  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)


  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:L2)
  U = TrialFESpace(V)

  T = TestFESpace(Ω, ReferenceFE(raviart_thomas,Float64,p); conformity=:Hdiv)
  S = TrialFESpace(T)

  Λ = ConstantFESpace(model)
  M = TrialFESpace(Λ)

  X = MultiFieldFESpace([S,U,M])
  Y = MultiFieldFESpace([T,V,Λ])

  ### force compatibility
  sigma_ex(x) = gradient(uex)(x)

  _X = MultiFieldFESpace([S,M])
  _Y = MultiFieldFESpace([T,Λ])

  biformS((s,μ),(t,λ)) = ∫( s⋅t )dΩ + ∫( divergence(s)*λ )dΩ  + ∫( divergence(t)*μ )dΩ
  liformS((t,λ)) = ∫( sigma_ex ⋅ t )dΩ
  op = AffineFEOperator(biformS,liformS,_X,_Y)
  sigma_exh,μh = solve(LUSolver(),op)

  # sum(∫((sigma_exh-sigma_ex)⊙(sigma_exh-sigma_ex))dΩ)


  ### # dual form -- with periodicity, force zero mean
  rhs = -1.0*divergence(sigma_exh)

  biformX((s,u,μ),(t,v,λ)) = (  ∫( s⋅t + divergence(t)*u )dΩ
                            + ∫( divergence(s)*v  )dΩ
                            + ∫(v*μ)dΩ + ∫(λ*u)dΩ
                      )
  liformY((t,v,λ)) = ∫( -(rhs*v) )dΩ  + ∫(λ*uex)dΩ

  op = AffineFEOperator(biformX,liformY,X,Y)
  sh,uh,μh = solve(LUSolver(),op)

  e = sum(∫((uh-uex)⊙(uh-uex))dΩ)
  println("Error = ", e)
  return e

end


p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns

domain = (0,1)

u(x) = 6*x[1]*(1-x[1]) # mean = 1.0 ∈ Ω = [0,1]
lab = "poly"

u(x) = sin(2*π*x[1]) + x[1]*(1-x[1])
lab = "trig"

# function u(x)
#   if x[1] < 0.5
#     return x[1]*(0.5-x[1]) + x[1]*(1-x[1])
#   else
#     return (x[1]-0.5)*(x[1]-1) + x[1]*(1-x[1])
#   end
# end
# lab = "pw"

errs = []
for n in collect(ns)
  partition = ntuple(x->n,Int(length(domain)/2))
  e = poisson_dual_form(domain,partition,p,degree,u)
  push!(errs,e)
end

plot(ns,errs,lw=3,markershape=:circle,label="Method 4: dual form")
plot!(yscale=:log10,xscale=:log10,framestyle=:box,legend=true,legendfontsize=11,show=true)
plot!(ns,1e-1dx.^4,lw=2,c=:black,ls=:dash,label="dx^4")
savefig(plotsdir()*"/poisson_dual_$(lab)_1D")



###############################################################################
########## low level

uex(x) = 6*x[1]*(1-x[1]) # mean = 1.0 ∈ Ω = [0,1]
uex(x) = sin(2*π*x[1]) + x[1]*(1-x[1])

p = 2
degree = 2*(p+1)
n = 16

model = CartesianDiscreteModel((0,1,0,1), (n,n),isperiodic=(true,true))

Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

writevtk(Ω,dir*"/poisson_dual",
        cellfields=["u"=>uex],append=false)

# dual form -- with periodicity, force zero mean
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:L2)
U = TrialFESpace(V)

T = TestFESpace(Ω, ReferenceFE(raviart_thomas,Float64,p); conformity=:Hdiv)
S = TrialFESpace(T)

Λ = ConstantFESpace(model)
M = TrialFESpace(Λ)

X = MultiFieldFESpace([S,U,M])
Y = MultiFieldFESpace([T,V,Λ])

### Sigma_exact
sigma_ex(x) = gradient(uex)(x)

_X = MultiFieldFESpace([S,M])
_Y = MultiFieldFESpace([T,Λ])

biformS((s,μ),(t,λ)) = ∫( s⋅t )dΩ + ∫( divergence(s)*λ )dΩ  + ∫( divergence(t)*μ )dΩ
liformS((t,λ)) = ∫( sigma_ex ⋅ t )dΩ
op = AffineFEOperator(biformS,liformS,_X,_Y)
sigma_exh,μh = solve(LUSolver(),op)

sum(∫((sigma_exh-sigma_ex)⊙(sigma_exh-sigma_ex))dΩ)


###
_rhs = -1.0*divergence(sigma_exh)
# rhs(x) = -1.0*divergence(sigma_ex)(x)
# rhs(x) = -1.0*laplacian(uex)(x)


biformX((s,u,μ),(t,v,λ)) = (  ∫( s⋅t + divergence(t)*u )dΩ
                            + ∫( divergence(s)*v  )dΩ
                            + ∫(v*μ)dΩ + ∫(λ*u)dΩ
                      )
liformY((t,v,λ)) = ∫( -(_rhs*v) )dΩ  + ∫(λ*uex)dΩ

op = AffineFEOperator(biformX,liformY,X,Y)
sh,uh,μh = solve(LUSolver(),op)

sum(∫((uh-uex)⊙(uh-uex))dΩ)

writevtk(Ω,dir*"/poisson_dual",
        cellfields=["u"=>uex, "uh"=>uh, "e"=>uex-uh,
        "s"=>sigma_ex,"sigma_h"=>sigma_exh ],append=false)
