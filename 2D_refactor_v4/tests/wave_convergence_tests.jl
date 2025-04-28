using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays


include("../src/initialise.jl")

#############
u_parametric(x) = zero(x)
p_parametric(x) = 1.0 + 0.001*x[1]

p = 2

######### cubed sphere model
coarse_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
model = Adaptivity.refine(coarse_model)
# ref_model = Adaptivity.refine(Adaptivity.refine(model))

manifold_model = model
panel_ids = get_panel_ids(manifold_model)

############# Surface metric and quadrature
order = 2*(p+1)
Ω_parametric = Triangulation(manifold_model)

m = Metric(cubedsphere,Ω_parametric)
dΩg = Measure(m,Ω_parametric,order)

############# FE problem
u_cf = CellField(u_parametric,Ω_parametric)
p_cf = CellField(p_parametric,Ω_parametric)


V = TestFESpace(coarse_model,ReferenceFE(raviart_thomas,Float64,order); conformity=:Hdiv)
Q = TestFESpace(manifold_model,ReferenceFE(lagrangian,Float64,order); conformity=:L2)

U = TrialFESpace(V)
P = TrialFESpace(Q)


f_u =  u_cf + surface_gradient(p_cf,m)
f_p = p_cf + surface_divergence(u_cf,m)

biform((u,p),(v,q)) =  ( ∫( u⋅v - surface_divergence(v,m)*p )dΩg
                       + ∫( p*q + surface_divergence(u,m)*q )dΩg )
liform((v,q))   =  ∫( f_u⋅v  + f_p⋅q )dΩg

op = AffineFEOperator(biform,liform,U,V)
ls = LUSolver()
uh = solve(ls,op)

A = get_matrix(op)
b = get_vector(op)


e_u = uh-u_cf
e_p = ph-p_cf
l2(e_u,dΩg)
l2(e_p,dΩg)

# plot in parametric space
writevtk(Ω_parametric,dir*"/possion_parametric",
        cellfields=["uh"=>uh,"u"=>u_cf, "e"=>e, "f"=>f_cf],append=false)


######## Plot in ambient space; map parametric -> ambient
ambient_model = AmbientDiscreteModel(manifold_model)
Ω_ambient = Triangulation(ambient_model)

cmap_parametric_2_ambient =  lazy_map(Reindex(Minv),panel_ids)

cf_mapped = lazy_map(Broadcasting(∘),get_data(u_cf),cmap_parametric_2_ambient)
u_ambient = CellData.similar_cell_field(u_cf,cf_mapped,Ω_ambient,PhysicalDomain() )

cf_mapped = lazy_map(Broadcasting(∘),get_data(uh),cmap_parametric_2_ambient)
uh_ambient = CellData.similar_cell_field(uh,cf_mapped,Ω_ambient,PhysicalDomain() )

e_ambient = u_ambient-uh_ambient
dΩ = Measure(Ω_ambient,order)
l2(e_ambient,dΩ)

writevtk(Ω_ambient,dir*"/possion_ambient",
        cellfields=["u_am"=>u_ambient, "uh_am"=>uh_ambient, "e_am"=>e_ambient],
        append=false)
