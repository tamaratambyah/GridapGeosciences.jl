ζ = 0.0
n_ref_lvls = 3
p_fe = 1
include("Williamson2Test.jl")

h = panel_to_cartesian(h₀(ζ))
vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
f = panel_to_cartesian(f₀(ζ))
b = panel_to_cartesian(_topography)

models  = get_refined_models(n_ref_lvls)
panel_model = models[1]


using GridapSolvers
ls = CGSolver(JacobiLinearSolver();rtol=1-12,verbose=1,name="diagnostic_solver")#

# get the ranks to help with storing/saving solution
ranks = get_ranks(panel_model)


## finite element solver
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

cell_geo_map = geo_map_func(Ω_panel)
dir = datadir("Distributed")
panel_cfs = [covarient_basis_cf⋅xh0[1], xh0[2]]
cellfields = map((x,y) -> x=>y, ["uh","ph"],panel_cfs)
writevtk(Ω_panel,dir*"/solT_0.vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)

cor_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
gravity = _g


# mectrics required in weak forms
metric_cf = CellField(analytic_metric,Ω_panel)
meas_cf = CellField(sqrtg,Ω_panel)
grad_meas_cf = CellField(grad_meas,Ω_panel)


#### DIAGNOSTIC VARIABLES


# vorticity
perp_matrix_cf = CellField(analytic_perp_matrix,Ω_panel)
aq(q,w) = ∫( q*xh0[2]*w*meas_cf  )dΩ
bq(w) = ∫( cor_cf*w*meas_cf  )dΩ + ∫( (perp_matrix_cf⋅xh0[1])⋅∇(w)  )dΩ
op = AffineFEOperator(aq,bq,H,R)
qh = solve(ls,op)

## mass flux
aF(F,v) =  ∫( (F⋅ (metric_cf⋅v))*meas_cf )dΩ
bF(v) = ∫( xh0[2]*(xh0[1]⋅(metric_cf⋅v))*meas_cf   )dΩ
op = AffineFEOperator(aF,bF,U,V)
Fh = solve(ls,op)

## bernouli poetnail
aΦ(Φ,ψ) = ∫( Φ*ψ*meas_cf  )dΩ
bΦ(ψ) = ∫( gravity*(xh0[2]+b_cf)*ψ*meas_cf  )dΩ + ∫( 0.5*( xh0[1] ⋅(metric_cf⋅xh0[1]) )ψ*meas_cf  )dΩ
op = AffineFEOperator(aΦ,bΦ,P,Q)
Φh = solve(ls,op)


A((q,F,Φ),(w,v,ψ)) = aq(q,w) + aF(F,v) + aΦ(Φ,ψ)
B((w,v,ψ)) = bq(w) + bF(v) + bΦ(ψ)
op = AffineFEOperator(A,B,X_diag,Y_diag)
yh = solve(ls,op)
_qh,_Fh,_Φh = yh

l2(qh-_qh,dΩ)
l2(Fh-_Fh,dΩ)
l2(Φh-_Φh,dΩ)

res((q,F,Φ),(w,v,ψ)) = A((q,F,Φ),(w,v,ψ)) - B((w,v,ψ))
jac((q,F,Φ),(dq,dF,dΦ),(w,v,ψ)) = A((dq,dF,dΦ),(w,v,ψ))
op = FEOperator(res,jac,X_diag,Y_diag)
nls = GridapSolvers.NonlinearSolvers.NewtonSolver(ls,verbose=true)
yh = solve(nls,op)
_qh,_Fh,_Φh = yh

l2(qh-_qh,dΩ)
l2(Fh-_Fh,dΩ)
l2(Φh-_Φh,dΩ)


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
nls = GridapSolvers.NonlinearSolvers.NewtonSolver(ls,verbose=true)
qh,Fh,Φh = solve(nls,_opFE)
vort = qh*xh0[2] - cor_cf


l2(qh-_qh,dΩ)
l2(Fh-_Fh,dΩ)
l2(Φh-_Φh,dΩ)


# equation for depth:
mass(t,(dut,dpt),(v,r)) = ∫( (dut⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( (dpt*r)*meas_cf )dΩ

res_p(((u,p),(q,F,Φ)),(v,r)) = ∫( r*(F⋅grad_meas_cf + meas_cf*(∇⋅F) )  )dΩ

res_u(((u,p),(q,F,Φ)),(v,r)) = (  ∫( q*( (perp_matrix_cf⋅F) ⋅(metric_cf ⋅v))   )dΩ
                              + ∫( -τ*(u⋅∇(q))*( (perp_matrix_cf⋅F) ⋅(metric_cf ⋅v))   )dΩ
                              - ∫( Φ*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
                  )

res_x(t,((u,p),(q,F,Φ)),(v,r)) = res_u(((u,p),(q,F,Φ)),(v,r)) + res_p(((u,p),(q,F,Φ)),(v,r))
jac_x(t,((u,p),(q,F,Φ)),(du,dp),(v,r)) =  ∫( -τ*(du⋅∇(q))*( (perp_matrix_cf⋅F) ⋅(metric_cf ⋅v))   )dΩ
jac_xt(t,((u,p),(q,F,Φ)),(dut,dpt),(v,r)) =  ∫( (dut⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( (dpt*r)*meas_cf )dΩ


opT = TransientSemilinearFEOperator(mass,res_x,(jac_x,jac_xt),X_prog,Y_prog)
opFE = FEOperator(res_y,jac_y,X_diag,Y_diag)
opDAE = DAEFEOperator(opT,opFE,ls)

CFL = 0.1
t0, tF = 0.0, _tF
_dt = dx(nc(panel_model))*CFL/(p_fe*sqrt(gravity*_H_0))
dt = floor(_dt, sigdigits=1)

τ = dt/2

ls_ode = CGSolver(JacobiLinearSolver();rtol=1-8,verbose=1,name="ode_solver")#
ls_ode.log.depth = 3
ode_solver = RungeKutta(ls_ode,ls_ode,dt,:EXRK_SSP_3_3)

solT  = solve(ode_solver,opDAE,t0,tF,xh0)
it = iterate(solT)


while !isnothing(it)
  data, state = it
  t, xh = data

  it = iterate(solT, state)
end
