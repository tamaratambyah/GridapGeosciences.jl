using Gridap
using GridapSolvers
using GridapGeosciences
using LinearAlgebra
using DrWatson

include("createpvd.jl")
p0(t) = x -> 1.0 + 0.0001*exp(-5*((1-x[1])^2+(-x[2])^2+(-x[3])^2))
u0(t) = x -> VectorValue(0.0 , 0.0, 0.0 )

n = 4
p = 0
degree = 6# 2*(p+1)

out_dir = datadir("wave_transient_n$(n)_p$(p)_highlevel")
!isdir(out_dir) && mkdir(out_dir)
pvd = createpvd(out_dir)


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

mass(t,(dtu,dtp),(v,q)) = ∫( dtu⋅v + dtp*q )dΩ
res(t,(u,p),(v,q)) =  -1.0*∫( (DIV(v))*p)dω + ∫( (DIV(u))*q )dω #∫( ∂t(u)⋅v + ∂t(p)*q )dΩ
jac(t,(u,p),(du,dp),(v,q)) = -1.0*∫( (DIV(v))*dp)dω + ∫( (DIV(du))*q )dω
jac_t(t,(u,p),(dtu,dtp),(v,q)) = ∫( dtu⋅v + dtp*q  )dΩ

opt = TransientSemilinearFEOperator(mass,res, (jac,jac_t), X, Y)
ls = LUSolver()
nls = NLSolver(ls, show_trace=true, method=:newton, iterations=10)


t0 = 0.0
tF = π
dt = tF/2000
CFL = dt/(1/n)

# solver = ThetaMethod(nls,dt,0.5)
solver = RungeKutta(ls,dt,:EXRK_SSP_2_2)
# solver = BackwardEuler(ls, dt)

xh0 = interpolate_everywhere([u0(0.0),p0(0.0)], X)
uh0,ph0 = xh0


pvd[0] = createvtk(Ω,out_dir*"/wave_equation_0.vtu",
                        cellfields=["u"=>uh0, "p"=>ph0],
                        append=false)


sol = solve(solver, opt, t0, tF, xh0)
it = iterate(sol)

@time while !isnothing(it)
  data, state = it
  t, xh = data
  uh, ph = xh
  println(t)
  pvd[t] = createvtk(Ω,out_dir*"/wave_equation_$t.vtu",
                      cellfields=["u"=>uh, "p"=>ph],
                                  append=false)
  it = iterate(sol, state)
end

make_pvd(out_dir,"wave_equation_","output",1)
