using Gridap
using Plots, LaTeXStrings
using DrWatson

include("../../src/initialise.jl")

# u(x) = VectorValue(cos(2*π*x[1]),cos(2*π*x[2]))
# h(x) = 1.0 + cos(2*π*x[1])*cos(2*π*x[2])



p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns


n = ns[2]
radius = 2

# Domain initialisation -- 1D
# domain = (0.0, 2*π)#, 0.0, 1.0)
# partition = (n, )# n)
# metric_func(x) = TensorValue{1}(radius^2)
# uex(x) = VectorValue(x[1]*(2π-x[1]))
# pex(x) = 1.0 + 0.01*x[1]*(2π-x[1])

# Domain initialisation -- 2D
domain = (-π/2,π/2,0,2π)
partition = (n,n)
metric_func(x) = TensorValue{2,2}(radius^2 ,0.0, 0.0, radius^2*(cos(x[1]))^2 )
uex(x) = VectorValue(0.0,x[2]*(2π-x[2]))
pex(x) = 1.0 +  0.01*(π/2+x[1])*(π/2-x[1])



##
model = CartesianDiscreteModel(domain, partition,isperiodic=ntuple(x->true,length(partition)))

Ω = Triangulation(model)
m = Metric(metric_func,Ω)

dΩ = Measure(Ω, degree)
dΩg =  Measure(m,Ω,degree)

ucf = CellField(uex,Ω)
pcf = CellField(pex,Ω)

writevtk(Ω,dir*"/wave",cellfields=["u"=>uex,"p"=>pex],append=false)


u0 = ucf + surface_gradient(pcf,m)
p0 = pcf + surface_divergence(ucf,m)

### FE problem
V = TestFESpace(model,ReferenceFE(raviart_thomas,Float64,p),conformity=:HDiv)
U = TrialFESpace(V)

Q = TestFESpace(model,ReferenceFE(lagrangian,Float64,p),conformity=:L2)
P = TrialFESpace(Q)

X = MultiFieldFESpace([U,P])
Y = MultiFieldFESpace([V,Q])

# Weak formulation for u
wave_biform((u,p),(v,q)) = ( ∫( u⋅v )dΩg
                   + ∫( -1.0*p*(wave_divergence(v,m))   )dΩ
                  + ∫( p*q )dΩg
                  + ∫( (surface_divergence(u,m))*q )dΩg
                  )
wave_liform((v,q)) = ∫( u0⋅v + p0*q  )dΩg

op = AffineFEOperator(wave_biform,wave_liform,X,Y)
uh, ph = solve(LUSolver(),op)

# Error
eu =  sum(∫((uh-uex)⊙(uh-uex))dΩ)
ep = sum(∫((ph-pex)⊙(ph-pex))dΩ)

eu_g =  sum(∫((uh-uex)⊙(uh-uex))dΩg)
ep_g = sum(∫((ph-pex)⊙(ph-pex))dΩg)

writevtk(Ω,dir*"/wave",
    cellfields=["u"=>uex,"p"=>pex,
                "uh"=>uh, "ph"=>ph,
                "eu"=>uh-uex,"ep"=>ph-pex],append=false)
