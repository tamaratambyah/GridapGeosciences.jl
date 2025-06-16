using Gridap, Gridap.Geometry
using Plots
using LinearAlgebra
using DrWatson
using Test
dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)



function poisson_periodic_1D(domain,partition,p,degree,uex;lagrange=false,uzeromean=false)
  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain,partition,isperiodic=ntuple(x->true,length(partition))))

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

    function poisson_liformY((v,λ))
      if uzeromean # force uex to have zeromean
        return ∫( rhs*v )dΩ  + ∫(λ*uex)dΩ
      else # only force u to have zero mean
        return  ∫( rhs*v )dΩ
      end
    end

    op = AffineFEOperator(poisson_biformX,poisson_liformY,X,Y)
    uh,μh = solve(BackslashSolver(),op)

    # b = get_matrix(op)
    # println("Compatibility: ", b)

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

domain = (0,1) #,0,1)

# polynomial
# u(x) = 6*x[1]*(1-x[1]) # mean = 1.0 ∈ Ω = [0,1]
# lab = "poly"

# trig
u(x) = sin(2*π*x[1]) + 1.0 # mean = 1.0 ∈ Ω = [0,1]
lab = "trig"

# piecewise
# function u(x)
#   if x[1] < 0.5
#     return x[1]*(0.5-x[1]) + 1.0
#   else
#     return (x[1]-0.5)*(x[1]-1) + 1.0
#   end
# end
# lab = "pw"

errs = []
errs_lagrange = []
errs_lagrange_zeromean = []
for n in collect(ns)
  partition = ntuple(x->n,Int(length(domain)/2))

  e = poisson_periodic_1D(domain,partition,p,degree,u)
  eg = poisson_periodic_1D(domain,partition,p,degree,u;lagrange=true)
  egz = poisson_periodic_1D(domain,partition,p,degree,u;lagrange=true,uzeromean=true)
  push!(errs,e)
  push!(errs_lagrange,eg)
  push!(errs_lagrange_zeromean,egz)
end

plot(ns,errs,lw=3,markershape=:circle,label="Method 1: H1(zeromean)")
plot!(ns,errs_lagrange,lw=3,markershape=:xcross,label="Method 2: ∫uₕ  = 0 ")
plot!(ns,errs_lagrange_zeromean,lw=3,markershape=:rect,label="Method 3: ∫uₕ-u = 0")
plot!(yscale=:log10,xscale=:log10,framestyle=:box,legend=true,legendfontsize=11,show=true)
plot!(ns,1e-1dx.^6,lw=2,c=:black,ls=:dash,label="dx^6")
# savefig(plotsdir()*"/poisson_periodic_$(lab)_1D")

## force zero mean
uzeromean(x) = u(x) - 1.0
errs = []
errs_lagrange = []
errs_lagrange_zeromean = []
for n in collect(ns)
  e = poisson_periodic_1D(domain,n,p,degree,uzeromean)
  eg = poisson_periodic_1D(domain,n,p,degree,uzeromean;lagrange=true)
  egz = poisson_periodic_1D(domain,n,p,degree,u;lagrange=true,uzeromean=true)
  push!(errs,e)
  push!(errs_lagrange,eg)
  push!(errs_lagrange_zeromean,egz)
end
plot(ns,errs,lw=3,markershape=:circle,label="Method 1: H1(zeromean)")
plot!(ns,errs_lagrange,lw=3,markershape=:xcross,label="Method 2: ∫uₕ  = 0 ")
plot!(ns,errs_lagrange_zeromean,lw=3,markershape=:rect,label="Method 3: ∫uₕ-u = 0")
plot!(yscale=:log10,xscale=:log10,framestyle=:box,legend=true,legendfontsize=11,show=true)
plot!(title="u₀ = u - 1")
# plot!(ns,dx.^2,lw=2,c=:black,ls=:dash,label="dx^2")
plot!(ns,1e-1dx.^6,lw=2,c=:black,ls=:dash,label="dx^6")
savefig(plotsdir()*"/poisson_periodic_$(lab)_force_1D")

###############################################################################
########## 2D tests
