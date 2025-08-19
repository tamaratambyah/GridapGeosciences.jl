
_sdiv(vec::Function,p) = αβ ->  sqrtg(αβ)*( vec(p)(αβ))
surfdiv(vec::Function,p::Int) = αβ -> 1/sqrtg(αβ) * ( divergence(_sdiv(vec,p))(αβ) )
surfdiv(vec::Function) = p -> surfdiv(vec,p)

f =  f_XYZ
vec_func(p) = αβ -> analytic_inv_metric(αβ) ⋅ gradient(f(p))(αβ)



panel_model = coarse_parametric_model()
panel_model = Gridap.Adaptivity.refine(panel_model)
panel_model = Gridap.Adaptivity.refine(panel_model)

p_fe = 1
Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,4*p_fe)
panel_ids = get_panel_ids(panel_model)

f_panel_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
sigma_cf = panelwise_cellfield(vec_func,Ω_panel,panel_ids)
sdiv_cf =  panelwise_cellfield(surfdiv(vec_func),Ω_panel,panel_ids)
slap_panel_cf =  panelwise_cellfield(surflap(f),Ω_panel,panel_ids)
covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

rhs_cf = f_panel_cf + sdiv_cf
# rhs_cf = f_panel_cf + divergence(sigma_cf)



l2(sdiv_cf-slap_panel_cf,dΩ)


V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
U = TrialFESpace(V)

T = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
S = TrialFESpace(T)

Y = MultiFieldFESpace([V, T])
X = MultiFieldFESpace([U, S])


inv_metric_cf = CellField(analytic_inv_metric,Ω_panel)
meas_cf = CellField(sqrtg,Ω_panel)

grad_meas(αβ) = gradient(sqrtg)(αβ)
div_inv_metric(αβ) = divergence(analytic_inv_metric)(αβ)

biform0((u,s),(v,t)) = ∫(  u*(meas_cf*(div_inv_metric⋅t + tr(analytic_inv_metric⋅gradient(t))  )
                            + (inv_metric_cf ⋅t)⋅grad_meas    )  )dΩ

biform1((u,s),(v,t)) = ∫( (s⋅t)*meas_cf )dΩ + biform0((u,s),(v,t)) #∫( (divergence(meas_cf*(inv_metric_cf⋅t)  ) )*u  )dΩ
biform2((u,s),(v,t)) = ∫( (u*v)*meas_cf )dΩ + ∫( v*(s⋅grad_meas + meas_cf*(∇⋅s) ) )dΩ

biformX((u,s),(v,t)) = biform1((u,s),(v,t)) + biform2((u,s),(v,t))
liformX((v,t)) = ∫( (rhs_cf*v)*meas_cf )dΩ  #+ ∫( (sigma_cf⋅t)*meas_cf )dΩ

op = AffineFEOperator(biformX,liformX,X,Y)
uh,sh = solve(LUSolver(),op)

e = f_panel_cf - uh
l2(e,dΩ)

e = sigma_cf - sh
l2(e,dΩ)


gradu = covarient_basis_cf ⋅sh

cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)

panel_cfs = [uh, sh, f_panel_cf,e, rhs_cf, gradu ]
labels = ["u","s","u_ex","e","rhs", "graduh"]

cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/ambient_model",cellfields=cellfields,append=false,geo_map=cell_geo_map)


_biform(u,t) = ∫( (divergence(meas_cf*(inv_metric_cf⋅t)  ) )*u  )dΩ
assemble_matrix(_biform,U,T)






biformX(s,t) = ∫( (s⋅t)*sqrtg )dΩ
liformX(t) = ∫( -1.0*(divergence(sqrtg*(analytic_inv_metric⋅t)  ) )*f_panel_cf  )dΩ
op = AffineFEOperator(biformX,liformX,S,T)
sh = solve(LUSolver(),op)
e = sigma_cf - sh
l2(e,dΩ)
