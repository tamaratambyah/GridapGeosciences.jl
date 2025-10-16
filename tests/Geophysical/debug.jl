using DrWatson
using Gridap, GridapGeosciences

include("Williamson2Test.jl")
include("../convergence_tools.jl")

ls_diag = LUSolver()
ls_ode = LUSolver()

n_ref_lvls = 4
p_fe = 1
ζ = 0.0
ls = LUSolver()
models  = get_refined_models(n_ref_lvls)

panel_model = models[1]

h = panel_to_cartesian(h₀(ζ))
vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
f = panel_to_cartesian(f₀(ζ))
η = panel_to_cartesian(η₀(ζ))
b = panel_to_cartesian(topography)

panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,2*(p_fe+1))

R = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1)
H = TransientTrialFESpace(R)

Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
P = TransientTrialFESpace(Q)

V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
U = TransientTrialFESpace(V)

X_prog = TransientMultiFieldFESpace([U,P]) # u, p
Y_prog = MultiFieldFESpace([V,Q]) # u, p

X_diag = TransientMultiFieldFESpace([H,U,P]) # q, F, Φ
Y_diag = MultiFieldFESpace([R,V,Q]) # q, F, Φ


## initial conditions
covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
u_contra_h = interpolate(u_contra_cf,U)
u_proj_h = covarient_basis_cf ⋅ u_contra_h

h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
b_cf = panelwise_cellfield(b,Ω_panel,panel_ids)
h_h = interpolate(h_cf-b_cf,P)

xh0 = interpolate_everywhere([u_contra_h,h_h],X_prog(0.0))


cor_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
gravity = _g


# mectrics required in weak forms
metric_cf = CellField(analytic_metric,Ω_panel)
meas_cf = CellField(sqrtg,Ω_panel)
grad_meas_cf = CellField(grad_meas,Ω_panel)


#### DIAGNOSTIC VARIABLES
# vorticity
perp_matrix_cf = CellField(analytic_perp_matrix,Ω_panel)
resq(((u,p),(q,F,Φ)),(w,v,ψ)) = ∫( q*p*w*meas_cf  )dΩ - ∫( cor_cf*w*meas_cf  )dΩ - ∫( (perp_matrix_cf⋅u)⋅∇(w)  )dΩ

# mass flux
resF(((u,p),(q,F,Φ)),(w,v,ψ)) = ∫( (F⋅ (metric_cf⋅v))*meas_cf )dΩ - ∫( p*(u⋅(metric_cf⋅v))*meas_cf   )dΩ

# Bernoulli potential
resΦ(((u,p),(q,F,Φ)),(w,v,ψ)) = ∫( Φ*ψ*meas_cf  )dΩ - ∫( gravity*(p+b_cf)*ψ*meas_cf  )dΩ - ∫( 0.5*( u ⋅(metric_cf⋅u) )ψ*meas_cf  )dΩ

res_y(t,((u,p),(q,F,Φ)),(w,v,ψ)) = resq(((u,p),(q,F,Φ)),(w,v,ψ)) + resF(((u,p),(q,F,Φ)),(w,v,ψ)) + resΦ(((u,p),(q,F,Φ)),(w,v,ψ))
jac_y(t,((u,p),(q,F,Φ)),(dq,dF,dΦ),(w,v,ψ)) = ∫( dq*p*w*meas_cf  )dΩ + ∫( (dF⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( dΦ*ψ*meas_cf  )dΩ

_res_y((q,F,Φ),(w,v,ψ))  = res_y(0.0,(xh0,(q,F,Φ)),(w,v,ψ))
_jac_y((q,F,Φ),(dq,dF,dΦ),(w,v,ψ)) = jac_y(0.0,(xh0,(q,F,Φ)),(dq,dF,dΦ),(w,v,ψ))
_opFE = FEOperator(_res_y,_jac_y,X_diag,Y_diag)
using GridapSolvers
ranks = get_ranks(panel_model)
nls = GridapSolvers.NonlinearSolvers.NewtonSolver(ls_diag,verbose=i_am_main(ranks))
qh,Fh,Φh = solve(nls,_opFE)
vort = qh*xh0[2] - cor_cf

#### PROGNOSTIC VARIABLES

# equation for depth:
mass(t,(dut,dpt),(v,r)) = ∫( (dut⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( (dpt*r)*meas_cf )dΩ

res_p(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) = ∫( r*(F⋅grad_meas_cf + meas_cf*(∇⋅F) )  )dΩ

res_u(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) = (  ∫( q*( (perp_matrix_cf⋅F) ⋅(metric_cf ⋅v))   )dΩ
                              + ∫( -τ*( (q-q0)/dt )*( (perp_matrix_cf⋅F) ⋅(metric_cf ⋅v))   )dΩ
                              + ∫( -τ*(u⋅∇(q))*( (perp_matrix_cf⋅F) ⋅(metric_cf ⋅v))   )dΩ
                              - ∫( Φ*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
                  )

res_x(t,((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) = res_u(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) + res_p(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0))
jac_x(t,((u,p),(q,F,Φ)),(du,dp),(v,r),(q0,F0,Φ0)) =  ∫( -τ*(du⋅∇(q))*( (perp_matrix_cf⋅F) ⋅(metric_cf ⋅v))   )dΩ
jac_xt(t,((u,p),(q,F,Φ)),(dut,dpt),(v,r),(q0,F0,Φ0)) =  ∫( (dut⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( (dpt*r)*meas_cf )dΩ


opT = TransientSemilinearFEOperator(mass,res_x,(jac_x,jac_xt),X_prog,Y_prog)
opFE = FEOperator(res_y,jac_y,X_diag,Y_diag)
opDAE = DAEFEOperator(opT,opFE,ls_diag)

CFL = 0.1
t0, tF = 0.0, _tF
_dt = dx(nc(panel_model))*CFL/(p_fe*sqrt(gravity*_H_0))
dt = floor(_dt, sigdigits=1)

τ = dt/2

ode_solver = RungeKutta(ls_ode,ls_ode,dt,:EXRK_SSP_3_3)

solT  = solve(ode_solver,opDAE,t0,tF,xh0)
it = iterate(solT)




Enstropys = Float64[]
Energys = Float64[]
Masss = Float64[]

ens0 = sum(∫( (qh*qh*xh0[2])*meas_cf  )dΩ)
energy0 = sum(∫( (0.5*xh0[2]*( xh0[1] ⋅(metric_cf⋅xh0[1])) + 0.5*gravity*xh0[2]*xh0[2] )*meas_cf )dΩ)
mass0 = sum( ∫( xh0[2]*meas_cf )dΩ  )
push!(Enstropys,ens0)
push!(Energys,energy0)
push!(Masss,mass0)

return_vtk = true
dir = datadir("test_supg")
!isdir(dir) && mkdir(dir)
cell_geo_map = geo_map_func(Ω_panel)
if return_vtk
  panel_cfs = [covarient_basis_cf⋅xh0[1], xh0[2],qh,Fh,Φh,vort,b_cf]
  cellfields = map((x,y) -> x=>y, ["uh","ph","qh","Fh","Phih","vort","bt"],panel_cfs)
  writevtk(Ω_panel,dir*"/solT_0.vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)
end


Es_u = Float64[]
Es_p = Float64[]
e_u,e_p = 0.0, 0.0
push!(Es_u,e_u)
push!(Es_p,e_p)

counter = 1
while !isnothing(it)
  data, state = it
  t, xh = data
  odeopcache = state[2][5][2]
  yh = odeopcache.diagnostics

  uh,ph = xh
  qh,Fh,Φh = yh

  vort = qh*ph - cor_cf
  println(t)

  uh_proj = covarient_basis_cf ⋅ uh
  e_u = l2( (u_proj_h - uh_proj)*meas_cf,dΩ)
  e_p = l2((h_cf - ph)*meas_cf,dΩ)
  push!(Es_u,e_u)
  push!(Es_p,e_p)

  # ens = sum(∫( (qh*qh*xh[2])*meas_cf  )dΩ)
  # energy = sum(∫( (0.5*xh[2]*( xh[1] ⋅(metric_cf⋅xh[1])) + 0.5*gravity*xh[2]*xh[2] )*meas_cf )dΩ)
  # _mass = sum( ∫( xh[2]*meas_cf )dΩ  )

  # push!(Enstropys,ens)
  # push!(Energys,energy)
  # push!(Masss,_mass)

  if return_vtk  && (mod(counter,10) == 0)
    panel_cfs = [covarient_basis_cf⋅uh, ph,qh,Fh,Φh,vort]
    cellfields = map((x,y) -> x=>y, ["uh","ph","qh","Fh","Phih","vort"],panel_cfs)
    writevtk(Ω_panel,dir*"/solT_$t.vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end
  counter = counter + 1
  it = iterate(solT, state)
end


make_pvd(dir,"solT",1)
x = rand(10)
_x = copy(x)
PartitionedArrays.consistent!(x) |> fetch
x == _x
