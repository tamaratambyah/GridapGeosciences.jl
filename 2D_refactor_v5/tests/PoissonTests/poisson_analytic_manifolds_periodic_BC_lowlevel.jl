using Gridap, Gridap.Geometry
using Plots
using LinearAlgebra
using DrWatson
using Test

include("../../src/initialise.jl")
# include("poisson_helpers.jl")



p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
n = ns[4]

r = 1.0


# Initialisation -- 1D
domain = (0,2π)
partition = (n,)
metric_func(x) = TensorValue{1}(r^2)
uex(x) = x[1]*(2π-x[1]) # mean = 1.0 ∈ Ω = [0,1]
# uex(x) = sin(x[1]) + x[1]*(2π-x[1])
# function uex(x)
#   if x[1] < π
#     return x[1]*(π-x[1]) # + x[1]*(1-x[1])
#   else
#     return (x[1]-π)*(x[1]-2*π) #+ x[1]*(1-x[1])
#   end
# end




#### Initialisation -- 2D -
#### Error blows up due to singularity at the pole. Therefore, only test in 1D
# domain = (-π/2,π/2, 0,2π)
# partition = (n,n)
# metric_func(x) = TensorValue{2,2}(r^2 ,0.0, 0.0, r^2*(cos(x[1]))^2 )
# function uex(x)
#   if x[2] < π
#     return x[2]*(π-x[2])
#   else
#     return (x[2]-π)*(x[2]-2*π)
#   end
# end

#######
model = CartesianDiscreteModel(domain, partition, isperiodic=ntuple(x->true,length(partition)))
Ω = Triangulation(model)
m = Metric(metric_func,Ω)

dΩ = Measure(Ω,degree)
dΩg =  Measure(m,Ω,degree)

# check zero mean
sum(∫(uex)dΩ  )

# check compatibility
ucf = CellField(uex,Ω)
sum(∫( surface_laplacian(ucf,m))dΩ  )

# writevtk(Ω,dir*"/poisson",cellfields=["u"=>uex],append=false)

# L2 = FESpace(model,ReferenceFE(lagrangian,Float64,1), conformity=:L2)
# biform(u,v) = ∫(u*v)dΩ
# liform(v) = ∫(v*ucf )dΩ
# op = AffineFEOperator(biform,liform,L2,L2)
# _uh = solve(LUSolver(),op)

################################################################################
#### Method 4
#### Mixed form -- with lagrange multiplers
################################################################################

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
sigma_ex(x) = surface_gradient(uex,m)(x)
# _sigma_ex = surface_gradient(_uh,m)
# writevtk(Ω,dir*"/poisson",cellfields=["s"=>sigma_ex],append=false)
# writevtk(cubedsphere,Ω,dir*"/poisson",cellfields=["s"=>sigma_ex,"_s"=>_sigma_ex],append=false)

_X = MultiFieldFESpace([S,M])
_Y = MultiFieldFESpace([T,Λ])

biformS((s,μ),(t,λ)) = ∫( s⋅t )dΩ + ∫( surface_divergence(s,m)*λ )dΩ  + ∫( surface_divergence(t,m)*μ )dΩ
liformS((t,λ)) = ∫( sigma_ex ⋅ t )dΩ
op = AffineFEOperator(biformS,liformS,_X,_Y)
sigma_exh,μh = solve(LUSolver(),op)

# biformS(s,t) = ∫( s⋅t )dΩ
# liformS(t) = ∫( _sigma_ex ⋅ t )dΩ
# op = AffineFEOperator(biformS,liformS,S,T)
# sigma_exh = solve(LUSolver(),op)

# writevtk(cubedsphere,Ω,dir*"/poisson",cellfields=["s"=>sigma_ex,"es"=>sigma_ex-sigma_exh],append=false)


### dual form
_rhs = -1.0*surface_divergence(sigma_exh,m)

biformX((s,u,μ),(t,v,λ)) = (  ∫( s⋅t + wave_divergence(t,m)*u )dΩ
                            + ∫( surface_divergence(s,m)*v  )dΩ
                            + ∫(v*μ)dΩ + ∫(λ*u)dΩ
                      )
liformY((t,v,λ)) = ∫( -(_rhs*v) )dΩ  + ∫(λ*uex)dΩ

op = AffineFEOperator(biformX,liformY,X,Y)
sh,uh,μh = solve(LUSolver(),op)

#### Compute errors
e = uh-uex
sum(∫(e⊙e)dΩ)
sum(∫(e⊙e)dΩg)

# writevtk(Ω,dir*"/poisson",
#         cellfields=["u"=>uex,"uh"=>uh,"e"=>e,"s_exh"=>sigma_exh,"s_ex"=>sigma_ex,
#         "es"=>sigma_ex-sigma_exh],append=false)
