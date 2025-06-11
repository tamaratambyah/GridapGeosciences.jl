using Gridap, Gridap.Geometry
using Plots
using LinearAlgebra
using DrWatson
using Test
dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)



function poisson_periodic_1D(domain,n,p,degree,uex;lagrange=false)
  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,(n,),isperiodic=(true,)))

  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  rhs(x) = -1.0*(laplacian(uex)(x))

  ### FE problem - multiifield
  if lagrange
    V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
    U = TrialFESpace(V)
    Λ = ConstantFESpace(model)
    M = TrialFESpace(Λ)
    X = MultiFieldFESpace([U,M])
    Y = MultiFieldFESpace([V,Λ])
    poisson_biformX((u,μ),(v,λ)) = ∫( gradient(u)⋅gradient(v)  )dΩ  + ∫(v*μ)dΩ + ∫(λ*u)dΩ
    poisson_liformY((v,λ)) = ∫( rhs*v )dΩ + ∫(λ*uex)dΩ
    op = AffineFEOperator(poisson_biformX,poisson_liformY,X,Y)
    uh,μh = solve(BackslashSolver(),op)
    return sum(∫((uh-uex)⊙(uh-uex))dΩ)
  end


  ### FE problem -- single field
  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1,constraint=:zeromean)
  U = TrialFESpace(V)
  poisson_biform(u,v) = ∫( gradient(u)⋅gradient(v) )dΩ
  poisson_liform(v) =  ∫( rhs*v )dΩ
  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
  uh = solve(BackslashSolver(),op)
  return sum(∫((uh-uex)⊙(uh-uex))dΩ)
end



p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns

domain = (0,1)

# polynomial
u(x) = 6*x[1]*(1-x[1]) # mean = 1.0 ∈ Ω = [0,1]
lab = "poly"

# trig
u(x) = sin(2*π*x[1]) + 1.0 # mean = 1.0 ∈ Ω = [0,1]
lab = "trig"

errs = []
errs_lagrange = []
for n in collect(ns)
  e = poisson_periodic_1D(domain,n,p,degree,u)
  eg = poisson_periodic_1D(domain,n,p,degree,u;lagrange=true)
  push!(errs,e)
  push!(errs_lagrange,eg)
end

plot(ns,errs,lw=3,label="constraint")
plot!(ns,errs_lagrange,lw=3,c=:red,label="lagrange")
plot!(yscale=:log10,xscale=:log10,framestyle=:box,legend=true,show=true)
savefig(plotsdir()*"/poisson_periodic_$(lab)_1D")

## force zero mean
uzeromean(x) = u(x) - 1.0
errs = []
errs_lagrange = []
for n in collect(ns)
  e = poisson_periodic_1D(domain,n,p,degree,uzeromean)
  eg = poisson_periodic_1D(domain,n,p,degree,uzeromean;lagrange=true)
  push!(errs,e)
  push!(errs_lagrange,eg)
end
plot(ns,errs,lw=3,label="constraint")
plot!(ns,errs_lagrange,lw=3,c=:red,label="lagrange")
plot!(yscale=:log10,xscale=:log10,framestyle=:box,legend=true,show=true)
savefig(plotsdir()*"/poisson_periodic_$(lab)_force_1D")

###############################################################################
##### low level -- debug
pts = get_cell_points(Ω)
cf_u = CellField(u,Ω)
cf_uzeromean = CellField(uzeromean,Ω)

# test periodicity
@test cf_u(pts)[1][1] == cf_u(pts)[end][2]
@test cf_uzeromean(pts)[1][1] == cf_uzeromean(pts)[end][2]
