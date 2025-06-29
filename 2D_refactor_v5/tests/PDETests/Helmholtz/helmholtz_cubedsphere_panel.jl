using Gridap, Gridap.Geometry
using Plots, LaTeXStrings
using LinearAlgebra
using DrWatson
using Test

include("../../../src/initialise.jl")
include("../pde_helpers.jl")


p = 2
degree = 20
ns = [2^i for i = 2:6]
dx =   ( sqrt.( 4*π*RADIUS^2 ./ (6*ns.^2) ) )


global RADIUS = 1.0*sqrt(3.0)


function helmholtz_cubedsphere_panel(n,p,degree,uex)

  model = CartesianDiscreteModel((-π/4,π/4, -π/4,π/4), (n,n))
  Ω = Triangulation(model)
  m = Metric(cubedsphere,Ω)

  dΩ = Measure(Ω,degree)
  dΩg =  Measure(m,Ω,degree)

  Γ = BoundaryTriangulation(model)
  dΓ = Measure(Γ,degree)
  dΓg = Measure(m,Γ,degree)
  n_Γ = get_normal_vector(Γ)

  ucf = CellField(uex,Ω)

  rhs = ucf + 1.0*(surface_laplacian(ucf,m))
  h = surface_gradient(ucf,m)⋅n_Γ

  # check compatibility
  compat = sum( ∫(rhs)dΩ   )
  println("Compatibility: ", compat)

  #### FE Problem -- no lagrange multiplers
  V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
  U = TrialFESpace(V)

  poisson_biform(u,v) = ∫(u*v)dΩg -  ∫( surface_gradient(u,m)⋅gradient(v)  )dΩg
  poisson_liform(v) = ∫(  rhs*v )dΩg - ∫( v*h )dΓg
  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

  uh = solve(LUSolver(),op)

  #### Compute errors
  e = l2(uh-ucf,dΩ)
  eg = l2(uh-ucf,dΩg)
  return e,eg

end



function u1(x)
  if x[2] < 0.0
    return -x[2]*(x[2] + π/4)
  else
    return x[2]*(x[2] - π/4)
  end

end
u2(x) = cos(4*x[1])*sin(4*x[2])
function u3(x)
  θϕ = GnomonicField()(x)
  return cos(θϕ[1])*sin(θϕ[2])
end
function u4(x)
  θϕ = GnomonicField()(x)
  X = SigmaField(RADIUS)(θϕ)
  return X[1]*X[2]*X[3]
end

uex_funcs = Dict{Symbol,Any}()
uex_funcs[:u1] = u1
uex_funcs[:u2] = u2
uex_funcs[:u3] = u3
# uex_funcs[:u4] = u4

errs = []
errs_g = []
for (key, val) in uex_funcs
  for n in ns
    e, eg = helmholtz_cubedsphere_panel(n,p,degree,val)
    push!(errs,e)
    push!(errs_g,eg)
  end
end

leginf = map(x->string(x),collect(keys(uex_funcs)))

plot()
plot_error(ns,errs;leginf=leginf,ls=fill(:solid,length(uex_funcs)))
plot_error(ns,errs_g;ls=fill(:dash,length(uex_funcs)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)",legend=:topright)
plot!(ns,1e-1dx.^6,lw=2,c=:black,label="dx^6")
savefig(plotsdir()*"/cubedsphere_panel_convergence")

# consider u4 separately as p=3 FE space
errs = []
errs_g = []
p = 3
degree = 2*(p+1)
for n in ns
  e, eg = helmholtz_cubedsphere_panel(n,p,degree,u4)
  push!(errs,e)
  push!(errs_g,eg)
end
plot()
plot_error(ns,errs;leginf=["u4"],ls=fill(:solid,length(1)))
plot_error(ns,errs_g;ls=fill(:dash,length(1)))
plot!(xscale=:log10,yscale=:log10,framestyle=:box,
xlabel="n cells",ylabel="L2(u - uh)",legend=:topright)
plot!(ns,1e-1dx.^6,lw=2,c=:black,label="dx^6")
savefig(plotsdir()*"/cubedsphere_panel_convergence_u4")

##############################################################################
##############l low level
n = 16


model = CartesianDiscreteModel((-π/4,π/4, -π/4,π/4), (n,n))
panel_ids = map(x->1,collect(1:num_cells(model)))

Ω = Triangulation(model)
m = Metric(cubedsphere,Ω)

dΩ = Measure(Ω, degree)
dΩg =  Measure(m,Ω,degree)

Γ = BoundaryTriangulation(model)
dΓ = Measure(Γ,degree)
dΓg = Measure(m,Γ,degree)
n_Γ = get_normal_vector(Γ)

writevtk(Γ,dir*"/poisson",append=false)

writevtk(Γ,dir*"/poisson",cellfields=["n"=>n_Γ],append=false)

uθϕ(θϕ) = cos(θϕ[1])*sin(θϕ[2])

function u_scalar_latlon2parametric(p::Int,uθϕ::Function)
  function _u(αβ)
    latlon_panel1 = GnomonicField()(αβ)
    uθϕ(latlon_panel1)
  end
end

# uex(x) = cos(4*x[1])*sin(4*x[2])
ucf = CellField(u4,Ω)

cell_field = map(p->GenericField(u_scalar_latlon2parametric(p,uθϕ)),panel_ids)
ucf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())
sum( ∫( ucf )dΩg), sum( ∫( ucf )dΩ)
sum(∫(surface_laplacian(ucf,m))dΩg  ), sum(∫(surface_laplacian(ucf,m))dΩ  )
sum( ∫( ucf )dΩg + ∫(surface_laplacian(ucf,m))dΩg  )
sum( ∫( ucf )dΩ + ∫(surface_laplacian(ucf,m))dΩ  )
writevtk(Ω,dir*"/poisson_og",cellfields=["u"=>ucf],append=false)

writevtk(Γ,dir*"/poisson",cellfields=["n"=>n_Γ,"uh"=>surface_gradient(ucf,m)⋅n_Γ],append=false)

h = surface_gradient(ucf,m)⋅n_Γ
rhs = ucf + 1.0*(surface_laplacian(ucf,m))
sum(∫( rhs)dΩg)

#### FE Problem -- no lagrange multiplers
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

poisson_biform(u,v) = ∫(u*v)dΩg -  ∫( surface_gradient(u,m)⋅gradient(v)  )dΩg
poisson_liform(v) = ∫(  rhs*v )dΩg - ∫( v*h )dΓg
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

uh = solve(LUSolver(),op)

A = get_matrix(op)
b = get_vector(op)
sum(b)


#### Compute errors
e = l2(uh-ucf,dΩ)
eg = l2(uh-ucf,dΩg)
writevtk(Ω,dir*"/poisson",cellfields=["u"=>ucf,"uh"=>uh,"e"=>ucf-uh],append=false)


################################################################################
#### Full cubed sphere
################################################################################
manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)
p = 2
degree = 10


Ω = Triangulation(manifold_model)
Γ = BoundaryTriangulation(manifold_model)
writevtk(Γ,dir*"/CS",append=false)


m = Metric(cubedsphere,Ω)

dΩ = Measure(Ω, degree)
dΩg =  Measure(m,Ω,degree)

# dΓ = Measure(Γ,degree)
# dΓg = Measure(m,Γ,degree)
# n_Γ = get_normal_vector(Γ)


function uex(x)
  if x[1] < 0.0
    return -x[1]*(x[1] + π/4)
  else
    return x[1]*(x[1] - π/4)
  end

end


ucf = CellField(uex,Ω)

# check zero mean and compatibility
sum(∫(ucf)dΩ  ), sum(∫( surface_laplacian(ucf,m))dΩ  )
sum(∫(ucf)dΩg  ), sum(∫( surface_laplacian(ucf,m))dΩg  )


rhs = -1.0*surface_laplacian(ucf,m)
# h = surface_gradient(ucf,m)⋅n_Γ


## zero mean in FE space
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

poisson_biform(u,v) = ∫(u*v)dΩg -  ∫( surface_gradient(u,m)⋅gradient(v)  )dΩg
poisson_liform(v) = ∫(  rhs*v )dΩg #- ∫( v*h )dΓg
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

uh = solve(LUSolver(),op)


e = l2(uh-ucf,dΩ)
eg = l2(uh-ucf,dΩg)
writevtk(Ω,dir*"/poisson",cellfields=["u"=>uex,"uh"=>uh,"e"=>uex-uh],append=false)
