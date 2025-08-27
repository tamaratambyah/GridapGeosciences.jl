
depth(XYZ) = 1.0 + 0.1*exp(-( XYZ[2]^2 + XYZ[3]^2 ) )
velocity(XYZ) = VectorValue(-XYZ[2],XYZ[1],0.0)
coriolis(XYZ) = 2.0

h = panel_to_cartesian(depth)
vecX = velocity
vX = panel_to_cartesian(tangent_vec(vecX))
f = panel_to_cartesian(coriolis)


ζ,u0,ω = π/2, 0.1, 1e-5
h = panel_to_latlon(hWilliamson(ζ,u0,ω))
vecX = vec_cartesian_to_latlon(vWilliamson(ζ,u0,ω))
vX = panel_to_cartesian(tangent_vec(vecX))
f = panel_to_latlon(fWilliamson(ζ,u0,ω))

p_fe = 1

## make models
panel_model = coarse_parametric_model()
panel_model = Gridap.Adaptivity.refine(panel_model)
panel_model = Gridap.Adaptivity.refine(panel_model)

panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,2*(p_fe+1))

Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
P = TrialFESpace(Q)

V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
U = TrialFESpace(V)


Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])


covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
u_proj_cf = panelwise_cellfield(projection_v(vX),Ω_panel,panel_ids)
cor_cf = panelwise_cellfield(f,Ω_panel,panel_ids)

u_contra = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
u_cf = covarient_basis_cf ⋅ u_contra

cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
panel_cfs = [h_cf,u_proj_cf,cor_cf , u_cf]
labels = ["p","u_proj","f", "ucf"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/ambient_model",cellfields=cellfields,append=false,geo_map=cell_geo_map)


u_perp_contra = panelwise_cellfield(contra_v_perp(vX),Ω_panel,panel_ids)
u_perp = covarient_basis_cf ⋅ u_perp_contra

sgrad_cf = panelwise_cellfield(sgrad(h),Ω_panel,panel_ids)
sdiv_cf =  panelwise_cellfield(surfdiv(contra_v(vX)),Ω_panel,panel_ids)

pinvJ_cf = panelwise_cellfield(forward_pinv_jacobian,Ω_panel,panel_ids)



# manufacture rhs functions
rhs_scalar = h_cf + sdiv_cf
rhs_vector = u_cf + cor_cf*u_perp + sgrad_cf
rhs_con_vector = pinvJ_cf ⋅ rhs_vector # exact contravariant component

geo_balance = cor_cf*u_perp + sgrad_cf
geo_balance_con = pinvJ_cf⋅ geo_balance
e_geo_balance = sum(∫( geo_balance_con )dΩ)
@check vector_length(e_geo_balance) < 1e-12

cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
panel_cfs = [rhs_con_vector,rhs_scalar,rhs_vector,geo_balance_con ]
labels = ["f1","f2", "rhs_vec","geo"]

cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/ambient_model",cellfields=cellfields,append=false,geo_map=cell_geo_map)

detg_cf = CellField(detg,Ω_panel)
metric_cf = CellField(analytic_metric,Ω_panel)
meas_cf = CellField(sqrtg,Ω_panel)
grad_meas_cf = CellField(grad_meas,Ω_panel)

function vecPerp(u)
  # u   = (u1, u2)
  # u^T = (-u2, u1)
  VectorValue(-u[2],u[1])
end

Aperp = [0 -1
        1 0]
Rperp = TensorValue(Aperp)
Rperp_cf = CellField(Rperp,Ω_panel)

biform1((u,p),(v,q)) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( ( cor_cf*( (Rperp_cf⋅ u)⋅v))*detg_cf )dΩ - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
biform2((u,p),(v,q)) = ∫( (p*q)*meas_cf )dΩ + ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) )  )dΩ

biformX((u,p),(v,q)) = biform1((u,p),(v,q)) + biform2((u,p),(v,q))
liformX((v,q)) = ∫( rhs_con_vector⋅(metric_cf⋅v)*meas_cf )dΩ + ∫( (rhs_scalar*q)*meas_cf )dΩ

op = AffineFEOperator(biformX,liformX,X,Y)
uh,ph = solve(LUSolver(),op)

uh_proj = covarient_basis_cf ⋅ uh

e_u = l2( (u_proj_cf - uh_proj)*meas_cf,dΩ) # error in physical velocity u
e_p = l2((h_cf - ph)*meas_cf,dΩ) # error in depth


lvl = nref(nc(panel_model))
cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
panel_cfs = [ph, uh_proj, uh_proj-u_proj_cf,ph-h_cf]
labels = ["p","u_proj","eu","ep"]

cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/ambient_model",cellfields=cellfields,append=false,geo_map=cell_geo_map)
