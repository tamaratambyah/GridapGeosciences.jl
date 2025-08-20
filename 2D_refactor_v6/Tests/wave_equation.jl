
_sdiv(vec::Function,p) = αβ ->  sqrtg(αβ)*( vec(p)(αβ))
surfdiv(vec::Function,p::Int) = αβ -> 1/sqrtg(αβ) * ( divergence(_sdiv(vec,p))(αβ) )
surfdiv(vec::Function) = p -> surfdiv(vec,p)

f =  f_sin
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

l2(sdiv_cf-slap_panel_cf,dΩ)


V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
U = TrialFESpace(V)

T = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
S = TrialFESpace(T)

Y = MultiFieldFESpace([V, T])
X = MultiFieldFESpace([U, S])

metric_cf = CellField(analytic_metric,Ω_panel)
inv_metric_cf = CellField(analytic_inv_metric,Ω_panel)
meas_cf = CellField(sqrtg,Ω_panel)
grad_meas_cf = CellField(grad_meas,Ω_panel)

biform1((u,s),(v,t)) = ∫( (s⋅ (metric_cf⋅t))*meas_cf )dΩ + ∫( u*(t⋅grad_meas_cf + meas_cf*(∇⋅t) ) )dΩ
biform2((u,s),(v,t)) = ∫( (u*v)*meas_cf )dΩ + ∫( v*(s⋅grad_meas_cf + meas_cf*(∇⋅s) ) )dΩ

biformX((u,s),(v,t)) = biform1((u,s),(v,t)) + biform2((u,s),(v,t))
liformX((v,t)) = ∫( (rhs_cf*v)*meas_cf )dΩ

op = AffineFEOperator(biformX,liformX,X,Y)
uh,sh = solve(LUSolver(),op)

e_u = f_panel_cf - uh
l2(e_u,dΩ)

e_s = sigma_cf - sh
l2(e_s,dΩ)


gradu = covarient_basis_cf ⋅sh

cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)

panel_cfs = [uh, sh, f_panel_cf,e_u, e_s,rhs_cf, gradu ]
labels = ["u","s","u_ex","e_u","e_s", "rhs", "graduh"]

cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/ambient_model",cellfields=cellfields,append=false,geo_map=cell_geo_map)
