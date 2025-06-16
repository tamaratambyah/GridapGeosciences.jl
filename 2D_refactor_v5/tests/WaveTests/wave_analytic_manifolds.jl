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
# pex(x) = 1.0 + 0.1*x[1]*(1-x[1]) + 0.1*x[2]*(1-x[2])

n = ns[2]
# Domain initialisation
domain = (0.0, 1.0)#, 0.0, 1.0)
partition = (n, )# n)

# metric_func(x) = TensorValue{1}(4)
metric_func(x) = TensorValue{1}(4*x[1]^2 + 1)
# metric_func(x) = TensorValue{1}(9*x[1]^4 + 12*x[1]^3 + 10*x[1]^2 + 4*x[1] + 2)

# metric_func(x) = TensorValue{2,2}(2,1,1,2 )
# metric_func(x) = TensorValue{2,2}(1+4*x[1]^2,4*x[1]*x[2],4*x[1]*x[2],1+4*x[2]^2 )

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
V = TestFESpace(model,ReferenceFE(raviart_thomas,Float64,p),conformity=:HDiv)
U = TrialFESpace(V)

W = TestFESpace(model,ReferenceFE(lagrangian,Float64,p),conformity=:L2)
R = TrialFESpace(W)

X = MultiFieldFESpace([U,R])
Y = MultiFieldFESpace([V,W])


function wave_divergence(v::CellField,m::Metric)
  # w =  m.inv_metric ⋅ v
  w = v⋅ m.inv_metric
  # gradient(w)*m.sq_meas + gradient(m.sq_meas) ⋅w
  m.sq_meas*divergence(w) +  gradient(m.sq_meas)⋅w


end

# Weak formulation for u
wave_biform((u,p),(v,q)) = ( ∫( u⋅v )dΩg
                   + ∫( -1.0*p*(wave_divergence(v,m))*1/m.sq_meas   )dΩg
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

tcf = Operation(transpose)(ucf)

pt = Point(0.5)
tcf(pt)⋅m.inv_metric(pt)
m.inv_metric(pt) ⋅ tcf(pt)

gradient(ucf ⋅ m.inv_metric)(pt)
gradient(w)*m.sq_meas + gradient(m.sq_meas) ⋅w


u(x) = VectorValue(x[1]*x[2],x[2])
# ∇u = gradient(u)(Point(2.0,1.0))
# ∇u[1,1] + ∇u[2,1]
# ∇u[1,2] + ∇u[2,2]


# prod =  (m.inv_metric ⋅ ucf)* m.sq_meas
# grad =  gradient(m.sq_meas)⋅  transpose(m.inv_metric ⋅ ucf)
# grad(Point(0.0))
# ( m.sq_meas *gradient(m.inv_metric ⋅ucf) )(Point(0.0))

# (m.inv_metric ⋅ ucf)(Point(0.0))

# gradient(m.sq_meas)(Point(0.0))
