using Gridap, Gridap.Geometry
using Plots
using LinearAlgebra
using DrWatson
using Test

include("../../../src/initialise.jl")



p = 2
degree = 20
ns = [2^i for i = 2:6]
n = 8


#######

global RADIUS = 1.0*sqrt(3.0)

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)
model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(model)


model = CartesianDiscreteModel((-π/4,π/4, -π/4,π/4), (n,n), isperiodic=(true,true))
panel_ids = map(x->1,collect(1:num_cells(model)))

Ω = Triangulation(model)
m = Metric(cubedsphere,Ω)

dΩ = Measure(Ω, degree)
dΩg =  Measure(m,Ω,degree)

uθϕ(θϕ) = cos(θϕ[1])*sin(θϕ[2])

function u_scalar_latlon2parametric(p::Int,uθϕ::Function)
  function _u(αβ)
    latlon_panel1 = GnomonicField()(αβ)
    uθϕ(latlon_panel1)
  end
end


cell_field = map(p->GenericField(u_scalar_latlon2parametric(p,uθϕ)),panel_ids)
ucf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())
sum( ∫( ucf )dΩg), sum( ∫( ucf )dΩ)
sum(∫(surface_laplacian(ucf,m))dΩg  ), sum(∫(surface_laplacian(ucf,m))dΩ  )
sum( ∫( ucf )dΩg + ∫(surface_laplacian(ucf,m))dΩg  )
sum( ∫( ucf )dΩ + ∫(surface_laplacian(ucf,m))dΩ  )
writevtk(Ω,dir*"/poisson",cellfields=["u"=>ucf],append=false)


rhs = ucf + 1.0*(surface_laplacian(ucf,m))
sum(∫( rhs)dΩg)

#### FE Problem -- no lagrange multiplers
V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

poisson_biform(u,v) = ∫(u*v)dΩg -  ∫( surface_gradient(u,m)⋅gradient(v)  )dΩg
poisson_liform(v) = ∫(  rhs*v )dΩg
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

uh = solve(LUSolver(),op)

A = get_matrix(op)
b = get_vector(op)
sum(b)


#### Compute errors
e = l2(uh-ucf,dΩ)
eg = l2(uh-ucf,dΩg)
writevtk(Ω,dir*"/poisson",cellfields=["u"=>ucf,"uh"=>uh,"e"=>ucf-uh,"_uh"=>uh+_ucf],append=false)
