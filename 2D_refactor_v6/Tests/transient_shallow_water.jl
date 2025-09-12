
panel_model = coarse_parametric_model()
panel_model = Gridap.Adaptivity.refine(panel_model)
panel_model = Gridap.Adaptivity.refine(panel_model)

p_fe = 1
CFL = 0.1
ζ = 0.0

h = panel_to_cartesian(h₀(ζ))
vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
f = panel_to_cartesian(f₀(ζ))
b = panel_to_cartesian(topography)


lvl = nref(nc(panel_model))
println("Refinement level: $lvl")

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
gravity = _g#1.0


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

#### PROGNOSTIC VARIABLES

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


ls = LUSolver()
opT = TransientSemilinearFEOperator(mass,res_x,(jac_x,jac_xt),X_prog,Y_prog)
opFE = FEOperator(res_y,jac_y,X_diag,Y_diag)
opDAE = DAEFEOperator(opT,opFE,ls)

t0, tF = 0.0, _tF#2*π
CFL = 0.1
_dt = dx(nc(panel_model))*CFL/(p_fe*sqrt(gravity*_H_0))
dt = 0.04

τ = dt/2

ode_solver = RungeKutta(ls,ls,dt,:EXRK_SSP_3_3)

solT  = solve(ode_solver,opDAE,t0,tF,xh0)

dir = datadir("Transient_shallow_water_W5_long")
!isdir(dir) && mkdir(dir)
labels = ["uh","ph","bt","h"]
panel_cfs = [covarient_basis_cf⋅xh0[1], xh0[2],b_cf,h_cf]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)

writevtk(Ω_panel,dir*"/solT_0.vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)

it = iterate(solT)

while !isnothing(it)
  data, state = it
  t, xh = data
  odeopcache = state[2][5][2]
  yh = odeopcache.diagnostics

  uh,ph = xh
  qh,Fh,Φh = yh

  vort = qh*ph - cor_cf
  println(t)

  panel_cfs = [covarient_basis_cf⋅uh, ph,qh,Fh,Φh,vort]
  cellfields = map((x,y) -> x=>y, ["uh","ph","qh","Fh","Phih","vort"],panel_cfs)
  writevtk(Ω_panel,dir*"/solT_$t.vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)

  it = iterate(solT, state)
end

make_pvd(dir,"solT",1)









uh_proj = covarient_basis_cf ⋅ uh
# e_u = l2( ( uh-u_contra_h  )*meas_cf,dΩ  )
e_u = l2( (u_proj_h - uh_proj)*meas_cf,dΩ) # error in physical velocity u
e_p = l2((h_cf - ph)*meas_cf,dΩ) # error in depth






if return_vtk
  lvl = nref(nc(panel_model))
  cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
  panel_cfs = [ph, h_cf,ph-h_cf,
              uh_proj, u_proj_h, uh_proj-u_proj_h,
               ]
  labels = ["ph","p","ep",
            "uh","u","eu",
              ]
  cellfields = map((x,y) -> x=>y, labels,panel_cfs)
  writevtk(Ω_panel,dir*"/ambient_model_nref$lvl",cellfields=cellfields,append=false,geo_map=cell_geo_map)
end

println(e_u, "; ", e_p, "; ",  e_η)





## Williamson2 convergence test
u0,ω, grav, H0 = 40/(6e6), 1e-5, 10, 3e3
u0,ω, grav, H0 = 0.1, 1e-5, 1.0, 1.0
n_ref_lvls = 4
williamson2_convergence_test(n_ref_lvls,true,u0,ω,grav,H0)
