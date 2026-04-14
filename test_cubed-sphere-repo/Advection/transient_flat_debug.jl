using DrWatson
using Test
using LinearAlgebra
using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.ODEs
using GridapGeosciences

dir = datadir("TransientTest")
!isdir(dir) && mkdir(dir)

# ut(t) = x -> x[1]^3 *t + x[2]^3*t
ut(t) = x -> x[1]^3 *t^3 + x[2]^3*t^3
u = TimeSpaceFunction(ut)
ft(t) = x -> ∂t(u)(t, x) - Δ(u)(t, x)
f = TimeSpaceFunction(ft)

# Geometry
domain = (0, 1, 0, 1)
partition = (5, 5)
model = CartesianDiscreteModel(domain, partition)

# FE spaces
p_fe = 3
reffe = ReferenceFE(lagrangian, Float64, p_fe)
V = FESpace(model, reffe, conformity=:H1, dirichlet_tags="boundary")
U = TransientTrialFESpace(V, u)

# Integration
Ω = Triangulation(model)
dΩ = Measure(Ω, 3*p_fe)

# FE operator
res(t, u, v) = ∫( ∂t(u)*v)dΩ +  ∫( ∇(u)⋅∇(v))dΩ - ∫( f(t)*v )dΩ
jac(t, u, du, v) = ∫( ∇(du)⋅∇(v))dΩ
jac_t(t, u, dut, v) =∫( dut*v)dΩ

opT = TransientFEOperator(res,(jac,jac_t),U,V)

# Initial conditions
t0 = 0.0
tF = 1.0
CFL = 0.1
dt = 0.001
# dt = 1/sqrt(num_cells(model))*CFL

uh0 = interpolate_everywhere(u(t0), U(t0))


ls = LUSolver()
nls = NLSolver(ls, show_trace=true, method=:newton, iterations=10)

tablus = [
  # :EXRK_SSP_3_3,
  :SDIRK_SSP_2_3,
  :SDIRK_3_3,
  :SDIRK_Norsett_3_4,
  :SDIRK_Crouzeix_3_4
]

dΩ_error = Measure(Ω, 9*p_fe)
Es = []
for tab in tablus
  println("$tab")
  solver = RungeKutta(nls, ls, dt, tab)
  solT = solve(solver, opT, t0, tF, uh0)

  errors = []
  for (t, uh) in solT
    e = u(t) - uh
    el2 = sqrt(sum(∫(e ⋅ e)dΩ_error))
    push!(errors,el2)
    # @test e_n < tol
    # writevtk(Ω,dir*"/solT_$t" * ".vtu", cellfields=["uh"=>uh,"u"=>u(t),"e"=>e],append=false)
  end
  push!(Es,errors)
  # make_pvd(dir,"solT",1)
end

using Plots
ts = collect(1:length(Es[1]))*dt
plot()
for (i,tab) in enumerate(tablus)
  plot!(ts,Es[i],lw=3,label="$tab")
end
plot!(xlabel="t",ylabel="L2(u(t)-uh(t))",
yscale=:log10,
legend=:bottomright)
plot!(show=true)
savefig(dir*"/cubic_big_dt")
