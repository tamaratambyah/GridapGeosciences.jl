using Gridap, Gridap.Geometry
using Plots
using LinearAlgebra
using DrWatson
using Test
dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns

domain = (0,1)
u(x) = 6*x[1]*(1-x[1]) # mean = 1.0 ∈ Ω = [0,1]
uzeromean(x) = u(x) - 1.0
laplacian_u = -12 # The laplacian is 12 for both u and uzeromean

uex(x) = uzeromean(x)

errs = []
for n in collect(ns)
  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,(n,),isperiodic=(true,)))

  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  pts = get_cell_points(Ω)
  cf_u = CellField(u,Ω)
  cf_uzeromean = CellField(uzeromean,Ω)

  # test periodicity
  @test cf_u(pts)[1][1] == cf_u(pts)[end][2]
  @test cf_uzeromean(pts)[1][1] == cf_uzeromean(pts)[end][2]

  # test laplacian
  rhs = CellField(-1.0*laplacian_u,Ω)
  _rhs(x) = -1.0*(laplacian(uex)(x))

  @test rhs(pts)[1][1] == _rhs(Point(domain[1]))
  @test rhs(pts)[end][2] == _rhs(Point(domain[2]))

  ### FE problem -- single field
  # V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1,constraint=:zeromean)
  # U = TrialFESpace(V)
  # poisson_biform(u,v) = ∫( gradient(u)⋅gradient(v) )dΩ
  # poisson_liform(v) =  ∫( rhs*v )dΩ
  # op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
  # uh = solve(LUSolver(),op)

  ### FE problem - multiifield
  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
  U = TrialFESpace(V)
  Λ = ConstantFESpace(model)
  M = TrialFESpace(Λ)
  X = MultiFieldFESpace([U,M])
  Y = MultiFieldFESpace([V,Λ])
  poisson_biform((u,μ),(v,λ)) = ∫( gradient(u)⋅gradient(v)  )dΩ  + ∫(v*μ)dΩ + ∫(λ*u)dΩ
  poisson_liform((v,λ)) = ∫( rhs*v )dΩ + ∫(λ*uex)dΩ
  op = AffineFEOperator(poisson_biform,poisson_liform,X,Y)
  uh,μh = solve(LUSolver(),op)


  e = sum(∫((uh-uex)⊙(uh-uex))dΩ)
  push!(errs,e)
  println("Errors: ", e)

end

plot(ns,errs,lw=3,label="singlefield")
plot!(ns,errs,lw=3,c=:red,label="multifield")
plot!(yscale=:log10,xscale=:log10,legend=true,show=true)
savefig(plotsdir()*"/poisson_periodic_1D")
