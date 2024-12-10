using Gridap
using GridapSolvers
using GridapGeosciences
using LinearAlgebra
using DrWatson

include("createpvd.jl")
p0(t) = x -> 1.0 + 0.0001*exp(-5*((1-x[1])^2+(-x[2])^2+(-x[3])^2))
u0(t) = x -> VectorValue(0.0 , 0.0, 0.0 )

n = 4
p = 1
degree = 4# 2*(p+1)

model = CubedSphereDiscreteModel(n)
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)
dω = Measure(Ω,degree,ReferenceDomain())
l2(e,dΩ) = sum(∫(e⊙e)dΩ)

rt_reffe = ReferenceFE(raviart_thomas, Float64, p)
lg_reffe = ReferenceFE(lagrangian, Float64, p)

V = FESpace(model, rt_reffe, conformity=:Hdiv)
U = TransientTrialFESpace(V)

Q = FESpace(model, lg_reffe; conformity=:L2)
P = TransientTrialFESpace(Q)

X = MultiFieldFESpace([U, P])
Y = MultiFieldFESpace([V, Q])


m(t, (dtu,dtp), (v,q)) = ∫( dtu⋅v + dtp*q )dΩ
a(t,(u, p), (v, q)) = ∫( - (DIV(v))*p)dω + ∫( (DIV(u))*q )dω
l(t,(v, q)) = ∫( VectorValue(0,0,0)⋅v + 0.0*q)dΩ

opt = TransientLinearFEOperator((m, a), l, X, Y, constant_forms=(true, true))


# PD = PatchDecomposition(model)
# P = GridapSolvers.VankaSolver(X,PD)

# P = GridapSolvers.VankaSolver(Y)
# ls = GMRESSolver(20;Pl=P,maxiter=1000,atol=1e-14,rtol=1.e-14,verbose=true)
ls = BackslashSolver()

dt = 0.05
solver = BackwardEuler(ls, dt)
t0 = 0.0
tF = 0.5

a0((u,p),(v,q)) = ∫( u⋅v + p*q )dΩ
l0((v,q)) = ∫( u0(t0)⋅v + p0(t0)*q )dΩ
op = AffineFEOperator(a0,l0,X,Y)
xh0 = solve(LUSolver(),op)



uh0,ph0 = xh0
pvd = createpvd(datadir("wave_equation"))
pvd[0] = createvtk(Ω,datadir("wave_equation")*"/wave_equation_0.vtu",
                        cellfields=["u"=>uh0, "p"=>ph0],
                        append=false)


sol = solve(solver, opt, t0, tF, xh0)
it = iterate(sol)

while !isnothing(it)
  data, state = it
  t, xh = data
  uh, ph = xh
  pvd[t] = createvtk(Ω,datadir("wave_equation")*"/wave_equation_$t.vtu",
                      cellfields=["u"=>uh, "p"=>ph],
                                  append=false)
  it = iterate(sol, state)
end

make_pvd(datadir("wave_equation"),"wave_equation",1)
