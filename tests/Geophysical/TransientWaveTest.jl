
using DrWatson
using Gridap
using GridapGeosciences
using Plots, LaTeXStrings

dir = datadir("Transient_wave_equation")
!isdir(dir) && mkdir(dir)
!isdir(plotsdir()) && mkdir(plotsdir())

panel_model = coarse_parametric_model()
panel_model = Gridap.Adaptivity.refine(panel_model)
panel_model = Gridap.Adaptivity.refine(panel_model)
panel_model = Gridap.Adaptivity.refine(panel_model)

## initial conditions
velocity(XYZ) = zero(XYZ)
depth(XYZ) = 1.0 + 0.01*exp(-5*((1-XYZ[1])^2+(0-XYZ[2])^2+(0-XYZ[3])^2))

## manipulate initial condition to sphere
h = panel_to_cartesian(depth)
vecX = velocity
vX = panel_to_cartesian(tangent_vec(vecX))

p_fe = 1

# transient parameters
t0, tF = 0.0, 2*π
CFL = 0.1
_dt = dx(nc(panel_model))*CFL/p_fe
dt = floor(_dt, sigdigits=1)

## finite element solver
panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,2*(p_fe+1))
cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)

Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
P = TransientTrialFESpace(Q)

V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
U = TransientTrialFESpace(V)

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

## initial conditions
h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
vec_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
xh0 = interpolate([vec_contra_cf,h_cf],X)

## transient weak form
metric_cf = CellField(analytic_metric,Ω_panel)
meas_cf = CellField(sqrtg,Ω_panel)
grad_meas_cf = CellField(grad_meas,Ω_panel)

mass(t, (dtu,dtp), (v,q)) = ∫( (v⋅ (metric_cf⋅ dtu) )*meas_cf )dΩ + ∫( (q*dtp)*meas_cf )dΩ
res(t,(u,p),(v,q)) =  ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) )  )dΩ - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
opT = TransientSemilinearFEOperator(mass, res, X, Y, constant_mass=true)

# solve with SSP RK 3
solver = RungeKutta(LUSolver(), LUSolver(), dt, :EXRK_SSP_3_3)
solT = solve(solver, opT, t0, tF, xh0)

## iterate solution
labels = ["uh","ph"]
panel_cfs = [xh0[1], xh0[2]]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/solT_0" * ".vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)

## casimirs to store
ts = Float64[]
ms = Float64[]
Es = Float64[]
s_divus = Float64[]
divus = Float64[]

_m0 = sum( ∫(  meas_cf*xh0[2] )dΩ)
_E0 = sum( ∫( 0.5*( xh0[1]⋅(metric_cf⋅ xh0[1]) + xh0[2]*xh0[2])*meas_cf )dΩ  )
_s_divu = sum(∫(  divergence(meas_cf*xh0[1]) )dΩ)
_divu = sum(∫(  divergence(xh0[1]) )dΩ)

push!(ts,0)
push!(ms,_m0)
push!(Es,_E0)
push!(s_divus,_s_divu)
push!(divus,_divu)

for (t, xh) in solT
  uh,ph = xh
  println(t)

  _m = sum( ∫(  meas_cf*ph )dΩ)
  _E = sum( ∫( 0.5*( uh⋅(metric_cf⋅ uh) + ph*ph)*meas_cf )dΩ  )
  _s_divu = sum(∫(   divergence(meas_cf*uh) )dΩ)
  _divu = sum(∫(  divergence(uh)  )dΩ)


  push!(ts,t)
  push!(ms,_m)
  push!(Es,_E)
  push!(s_divus,_s_divu)
  push!(divus,_divu)

  panel_cfs = [uh, ph]
  cellfields = map((x,y) -> x=>y, labels,panel_cfs)
  writevtk(Ω_panel,dir*"/solT_$t" * ".vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)
end


make_pvd(dir,"solT",1)

# plot relative error in mass and energy
ms_rel = abs.(ms.-ms[1])./ms[1]
Es_rel = abs.(Es.-Es[1])./Es[1]

plot()
plot!(ts[2:end],ms_rel[2:end],lw=3,label="mass")
plot!(ts[2:end],Es_rel[2:end],lw=3,label="energy")
plot!(yaxis=:log,xlabel="t",ylabel=L"|x_t-x_0|/x_0")
savefig(plotsdir()*"/wave_transient_conservation")

# plot surface divergence
plot()
plot!(ts[2:end],abs.(s_divus[2:end]),lw=3)
plot!(yaxis=:log,xlabel="t",ylabel="|s_div (u)|")
savefig(plotsdir()*"/wave_transient_s_div")

# plot panel divergence
plot()
plot!(ts[2:end],abs.(divus[2:end]),lw=3)
plot!(xlabel="t",ylabel="|div (u)|")
savefig(plotsdir()*"/wave_transient_panel_div")
