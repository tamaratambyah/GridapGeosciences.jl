using Gridap, Gridap.Geometry
using Plots, LaTeXStrings
using LinearAlgebra
using DrWatson

include("../../../src/initialise.jl")
include("../pde_helpers.jl")



p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns

global RADIUS = 1.0*sqrt(3.0)

uex(x) = VectorValue(0.0,(π/4+x[1])*(π/4-x[1]))
pex(x) = 1.0 +  0.01*(π/4+x[1])*(π/4-x[1])





###############################################################################
global RADIUS = 1.0*sqrt(3.0)

n = 16
p = 2
degree = 2*(p+1)

uex(x) = VectorValue(0.0,(π/4+x[1])*(π/4-x[1]))
pex(x) = 1.0 +  0.01*(π/4+x[1])*(π/4-x[1])

# uex(x) = VectorValue(0.0,cos(4*x[2]))
# pex(x) = 1.0 +  0.1*cos(4*x[1])


##
model = CartesianDiscreteModel((-π/4,π/4, -π/4,π/4),(n,n) )

Ω = Triangulation(model)
m = Metric(cubedsphere,Ω)


dΩ = Measure(Ω, degree)
dΩg =  Measure(m,Ω,degree)

Γ = BoundaryTriangulation(model)
m2 = Metric(cubedsphere,Γ)


dΓ = Measure(Γ,degree)
dΓg = Measure(m,Γ,degree)
n_Γ = get_normal_vector(Γ)

ucf = CellField(uex,Ω)
pcf = CellField(pex,Ω)

writevtk(Ω,dir*"/wave_cubed_sphere",cellfields=["u"=>ucf,"p"=>pcf],append=false)


u0 = ucf + surface_gradient(pcf,m)
p0 = pcf + surface_divergence(ucf,m)

### FE problem
V = TestFESpace(model,ReferenceFE(raviart_thomas,Float64,p),conformity=:HDiv)
U = TrialFESpace(V)

Q = TestFESpace(model,ReferenceFE(lagrangian,Float64,p),conformity=:L2)
P = TrialFESpace(Q)

X = MultiFieldFESpace([U,P])
Y = MultiFieldFESpace([V,Q])

# p_trial = get_trial_fe_basis(P)
# v_test = get_fe_basis(V)
# liform(v) = ∫( pcf*(m.inv_metric⋅v)⋅n_Γ )dΓg
# liform(v_test)


# Weak formulation for u
wave_biform((u,p),(v,q)) = ( ∫( u⋅v )dΩg
                   + ∫( -1.0*p*(wave_divergence(v,m))   )dΩ
                  + ∫( p*q )dΩg
                  + ∫( (surface_divergence(u,m))*q )dΩg
                  )
wave_liform((v,q)) = ∫( u0⋅v + p0*q  )dΩg -  ∫( pcf*(m2.inv_metric⋅v)⋅n_Γ )dΓg

op = AffineFEOperator(wave_biform,wave_liform,X,Y)
uh, ph = solve(LUSolver(),op)

# Error
eu =  sum(∫((uh-ucf)⊙(uh-ucf))dΩ)
ep = sum(∫((ph-pcf)⊙(ph-pcf))dΩ)

eu_g =  sum(∫((uh-ucf)⊙(uh-ucf))dΩg)
ep_g = sum(∫((ph-pcf)⊙(ph-pcf))dΩg)

writevtk(Ω,dir*"/wave_cubed_sphere",
    cellfields=["u"=>ucf,"p"=>pcf,
                "uh"=>uh, "ph"=>ph,
                "eu"=>uh-ucf,"ep"=>ph-pcf],append=false)


# split: u equation only
wave_biform(u,v) = ( ∫( u⋅v )dΩg
                  )
wave_liform(v) = ∫( u0⋅v )dΩg + ∫( pcf*(wave_divergence(v,m))   )dΩ -  ∫( pcf*(m2.inv_metric⋅v)⋅n_Γ )dΓg

op = AffineFEOperator(wave_biform,wave_liform,U,V)
uh = solve(LUSolver(),op)
eu =  sum(∫((uh-ucf)⊙(uh-ucf))dΩ)
eu_g =  sum(∫((uh-ucf)⊙(uh-ucf))dΩg)

# split: p equation only
wave_biform(p,q) = (  ∫( p*q )dΩg
                  )
wave_liform(q) = ∫( p0*q  )dΩg  - ∫( (surface_divergence(ucf,m))*q )dΩg
op = AffineFEOperator(wave_biform,wave_liform,P,Q)
ph = solve(LUSolver(),op)

ep = sum(∫((ph-pcf)⊙(ph-pcf))dΩ)
ep_g = sum(∫((ph-pcf)⊙(ph-pcf))dΩg)
