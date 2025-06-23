using Gridap, Gridap.Geometry
using Plots, LaTeXStrings
using LinearAlgebra
using DrWatson
using Test

include("../../../src/initialise.jl")
include("../poisson_helpers.jl")


p = 2
degree = 20
ns = [2^i for i = 2:6]
dx = 1 ./ ns

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


uex_funcs = Dict{Symbol,Any}()
uex_funcs[:u1] = u1
uex_funcs[:u2] = u2
uex_funcs[:u3] = u3

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
# ucf = CellField(uex,Ω)

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
