using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays


include("../src/initialise.jl")

######## mappings
R = lazy_map(x-> PanelRotationField(r1p[x]), 1:6)
_R = lazy_map(x-> PanelRotationField(rp1[x]), 1:6)
γ = GnomonicField()
γinv = InvGnomonicField()
σ = SigmaField(r)

M = lazy_map(x->    R[x]∘ σ ∘ γ, 1:6)
Minv = lazy_map(x-> γinv ∘ σ ∘ _R[x], 1:6)

#########

coarse_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
model = Adaptivity.refine(coarse_model)
ref_model = Adaptivity.refine(model)

p = 2
order = 2*p+1

u_ambient(x) = x[1]*x[2]*x[3]

manifold_model = ref_model
ambient_model = AmbientDiscreteModel(manifold_model)
panel_ids = get_panel_ids(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)

pts_ambient = get_cell_points(Ω_ambient)
pts_parametric = get_cell_points(Ω_parametric)

######## ambient -> parametric
cf_ambient  = CellField(u_ambient,Ω_ambient)
cell_field = get_data(cf_ambient)
cell_invmap =  lazy_map(Reindex(M),panel_ids)

cell_field_phys = lazy_map(Broadcasting(∘),cell_field,cell_invmap)
cf_parametric = CellData.similar_cell_field(cf_ambient,cell_field_phys,Ω_parametric,PhysicalDomain() )
cf_parametric(pts_parametric)

writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)
writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)


m = Metric(cubedsphere,Ω_parametric)
dΩg = Measure(m,Ω_parametric,order)
l2(e,dΩg) = sum(∫(e⊙e)dΩg)

### not working
f_cf = -1.0*surface_laplacian(cf_parametric,m)



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
