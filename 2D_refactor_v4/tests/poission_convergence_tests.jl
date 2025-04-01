using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../src/initialise.jl")

manifold_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
model = Adaptivity.refine(manifold_model)
p = 1
order = 2*p+1

u(x) = x[1] #*(a-x[1]) + x[2]*(a-x[2])

Ω = Triangulation(model)
m = Metric(cubedsphere,Ω)
dΩg = Measure(m,Ω,order)
l2(e,dΩg) = sum(∫(e⊙e)dΩg)


u_cf = CellField(u,Ω)
f_cf = -1.0*surface_laplacian(u_cf,m)

V = FESpace(model,ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

bilinear(u,v) = ∫(  surface_gradient(v,m)⋅ surface_gradient(u,m)  )dΩg
linear(v)   = ∫( v*f_cf )dΩg

op = AffineFEOperator(bilinear,linear,U,V)
ls = LUSolver()
uh = solve(ls,op)

e = uh-u_cf

dΩ = Measure(Ω,order)

l2(e,dΩ)

writevtk(Ω,dir*"/possion",cellfields=["uh"=>uh,"u"=>u_cf, "e"=>e],append=false)
# writevtk(get_ambient_grid(get_grid(model)),dir*"/possion",cellfields=["uh"=>uh,"u"=>u_cf, "e"=>e],append=false)
