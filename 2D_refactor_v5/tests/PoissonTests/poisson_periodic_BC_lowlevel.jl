using Gridap, Gridap.Geometry
using Plots
using LinearAlgebra
using DrWatson
using Test
dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


p = 2
degree = 2*(p+1)

# uex(x) = 6*x[1]*(1-x[1]) # mean = 1.0 ∈ Ω = [0,1]
# uex(x) = sin(2*π*x[1]) + x[1]*(1-x[1])
function uex(x)
  if x[1] < 0.5
    return x[1]*(0.5-x[1]) # + x[1]*(1-x[1])
  else
    return (x[1]-0.5)*(x[1]-1) #+ x[1]*(1-x[1])
  end
end


model = CartesianDiscreteModel((0,1,0,1), (8,8), isperiodic=(true,true))
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

# check zero mean
sum(∫(uex)dΩ  )

# check compatibility
sum(∫( laplacian(uex))dΩ  )

writevtk(Ω,dir*"/poisson",cellfields=["u"=>uex],append=false)

################################################################################
#### Method 1
#### primial form -- no lagrange multiplers
################################################################################
rhs(x) = -1.0*(laplacian(uex)(x))

#### FE Problem -- no lagrange multiplers
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1, constraint=:zeromean)
U = TrialFESpace(V)

poisson_biform(u,v) = ∫( gradient(u)⋅gradient(v)  )dΩ
poisson_liform(v) = ∫(  rhs*v )dΩ
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

A = get_matrix(op)
b = get_vector(op)

# evals = eigvals(Array(A))
# println(A*(3*ones(size(b))))
sum(b)

uh = solve(LUSolver(),op)
sum( ∫(uh )dΩ  )

#### Compute errors
e = uh-uex
sum(∫(e⊙e)dΩ)

writevtk(Ω,dir*"/poisson",
        cellfields=["u"=>uex,"uh"=>uh,"e"=>e],append=false)


################################################################################
#### Method 2 / 3
#### primial form -- with lagrange multiplers
################################################################################
rhs(x) = -1.0*(laplacian(uex)(x))

# dual form -- with periodicity, force zero mean
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

Λ = ConstantFESpace(model)
M = TrialFESpace(Λ)

X = MultiFieldFESpace([U,M])
Y = MultiFieldFESpace([V,Λ])

###
poisson_biformX((u,μ),(v,λ)) = ∫( gradient(u)⋅gradient(v)  )dΩ  + ∫(v*μ)dΩ + ∫(λ*u)dΩ
poisson_liformY((v,λ)) =  ∫( rhs*v )dΩ  + ∫(λ*uex)dΩ

op = AffineFEOperator(poisson_biformX,poisson_liformY,X,Y)
uh,μh = solve(LUSolver(),op)

#### Compute errors
e = uh-uex
sum(∫(e⊙e)dΩ)

writevtk(Ω,dir*"/poisson",
        cellfields=["u"=>uex,"uh"=>uh,"e"=>e],append=false)



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
sigma_ex(x) = gradient(uex)(x)

_X = MultiFieldFESpace([S,M])
_Y = MultiFieldFESpace([T,Λ])

biformS((s,μ),(t,λ)) = ∫( s⋅t )dΩ + ∫( divergence(s)*λ )dΩ  + ∫( divergence(t)*μ )dΩ
liformS((t,λ)) = ∫( sigma_ex ⋅ t )dΩ
op = AffineFEOperator(biformS,liformS,_X,_Y)
sigma_exh,μh = solve(LUSolver(),op)

# sum(∫((sigma_exh-sigma_ex)⊙(sigma_exh-sigma_ex))dΩ)


### dual form
_rhs = -1.0*divergence(sigma_exh)

biformX((s,u,μ),(t,v,λ)) = (  ∫( s⋅t + divergence(t)*u )dΩ
                            + ∫( divergence(s)*v  )dΩ
                            + ∫(v*μ)dΩ + ∫(λ*u)dΩ
                      )
liformY((t,v,λ)) = ∫( -(_rhs*v) )dΩ  + ∫(λ*uex)dΩ

op = AffineFEOperator(biformX,liformY,X,Y)
sh,uh,μh = solve(LUSolver(),op)

#### Compute errors
e = uh-uex
sum(∫(e⊙e)dΩ)

writevtk(Ω,dir*"/poisson",
        cellfields=["u"=>uex,"uh"=>uh,"e"=>e],append=false)
