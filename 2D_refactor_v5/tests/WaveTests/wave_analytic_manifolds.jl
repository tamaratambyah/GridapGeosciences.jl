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

uex(x) = VectorValue(x[1]*(1-x[1]))
pex(x) = 1.0 + 0.1*x[1]*(1-x[1])


# uex(x) = VectorValue(x[1]*(1-x[1]),x[2]*(1-x[2]))
# pex(x) = 1.0 + x[1]*(1-x[1]) + x[2]*(1-x[2])

n = ns[2]
# Domain initialisation
domain = (0.0, 1.0)#, 0.0, 1.0)
partition = (n, )#n)
metric_func(x) = TensorValue{1}(4)

model = CartesianDiscreteModel(domain, partition,isperiodic=(true,))

Ω = Triangulation(model)
m = Metric(metric_func,Ω)

dΩ = Measure(Ω, degree)
dΩg =  Measure(m,Ω,degree)

ucf = CellField(uex,Ω)
pcf = CellField(pex,Ω)

u0 = ucf + surface_gradient(pcf,m)
p0 = pcf + surface_divergence(ucf,m)

### FE problem
V = TestFESpace(model,
    ReferenceFE(raviart_thomas,Float64,p),
    conformity=:HDiv)
U = TrialFESpace(V)

W = TestFESpace(model,
    ReferenceFE(lagrangian,Float64,p),
    conformity=:L2)
R = TrialFESpace(W)

X = MultiFieldFESpace([U,R])
Y = MultiFieldFESpace([V,W])


function wave_divergence(v::CellField,m::Metric)
  _v =  m.inv_metric ⋅ v
  m.sq_meas*divergence(_v) +  gradient(m.sq_meas)⋅_v
end

# Weak formulation for u
wave_biform((u,p),(v,q)) = ( ∫( u⋅v )dΩg
                           + ∫( -1.0*p*(wave_divergence(v,m))   )dΩ
                  + ∫( p*q )dΩg
                  + ∫( (surface_divergence(u,m))*q )dΩg
                  )
wave_liform((v,q)) = ∫( u0⋅v + p0*q  )dΩ

op = AffineFEOperator(wave_biform,wave_liform,X,Y)
uh, ph = solve(op)

# Error
eu =  sum(∫((uh-uex)⊙(uh-uex))dΩ)
ep = sum(∫((ph-pex)⊙(ph-pex))dΩ)

eu_g =  sum(∫((uh-uex)⊙(uh-uex))dΩg)
ep_g = sum(∫((ph-pex)⊙(ph-pex))dΩg)
