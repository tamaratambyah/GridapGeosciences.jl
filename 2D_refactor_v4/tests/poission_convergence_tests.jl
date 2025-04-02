using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../src/initialise.jl")

coarse_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
model = Adaptivity.refine(coarse_model)
ref_model = Adaptivity.refine(model)

p = 2
order = 2*p+1

u_ambient(x) = x[1]*x[2]*x[3]

manifold_model = ref_model
ambient_model = AmbientDiscreteModel(ref_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)


pts_ambient = get_cell_points(Ω_ambient)
pts_parametric = get_cell_points(Ω_parametric)

######## ambient -> parametric
cf_ambient = change_domain(CellField(u_ambient,Ω_ambient),ReferenceDomain())
cf_parametric = GenericCellField(CellData.get_data(cf_ambient),Ω_parametric,ReferenceDomain())
_cf_parametric = change_domain(cf_parametric,PhysicalDomain())

writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric, "v"=>_cf_parametric],append=false)

m = Metric(cubedsphere,Ω_parametric)
dΩg = Measure(m,Ω_parametric,order)
l2(e,dΩg) = sum(∫(e⊙e)dΩg)

f_cf = -1.0*surface_laplacian(_cf_parametric,m)

### not working
grad = gradient(_cf_parametric)
evaluate(grad,pts_parametric.cell_phys_point[1])
grad(pts_parametric)



V = FESpace(manifold_model,ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

bilinear(u,v) = ∫(  surface_gradient(v,m)⋅ surface_gradient(u,m)  )dΩg
linear(v)   = ∫( v*f_cf )dΩg

op = AffineFEOperator(bilinear,linear,U,V)
ls = LUSolver()
uh = solve(ls,op)

e = uh-u_cf
l2(e,dΩg)

writevtk(Ω,dir*"/possion",cellfields=["uh"=>uh,"u"=>u_cf, "e"=>e],append=false)
# writevtk(get_ambient_grid(get_grid(model)),dir*"/possion",cellfields=["uh"=>uh,"u"=>u_cf, "e"=>e],append=false)
