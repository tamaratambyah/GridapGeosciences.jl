## make models
cube_model = coarse_cube_model(π/4,6)
cube_model = Gridap.Adaptivity.refine(cube_model)
cube_model = Gridap.Adaptivity.refine(cube_model)

panel_model = parametric_model(cube_model)
sphere_model = ambient_model(panel_model)

# f = f_sin
# f = panel_to_cartesian(fX)
f = panel_to_latlon(fWilliamson(π/2))

################################################################################
#### On panel
################################################################################
### Panel cell fields
Ω_panel = Triangulation(panel_model)
panel_ids = get_panel_ids(panel_model)

f_panel_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
contravarient_basis_cf = panelwise_cellfield(contravariant_basis,Ω_panel,panel_ids)
inv_metric_cf = CellField(analytic_inv_metric,Ω_panel)
slap_panel_cf =  panelwise_cellfield(surflap(f),Ω_panel,panel_ids)


### analytic gradient -- need to define new function to trigger automatric differentiation
gradf(p) = αβ -> gradient(f(p))(αβ)
gradf_cf = panelwise_cellfield(gradf,Ω_panel,panel_ids)

grad_covarient =  contravarient_basis_cf ⋅ gradf_cf
grad_contravarient = covarient_basis_cf ⋅  (inv_metric_cf ⋅ gradf_cf)

### interpolate f into FE space, and then differentiate
p_fe = 2
V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
U = TrialFESpace(V)

f_uh = interpolate(f_panel_cf,U)

grad_covarient_uh =  contravarient_basis_cf ⋅ gradient(f_uh)
grad_contravarient_uh = covarient_basis_cf ⋅  (inv_metric_cf ⋅ gradient(f_uh))

e_u = f_panel_cf - f_uh
e_con = grad_contravarient-grad_contravarient_uh
e_cov = grad_covarient-grad_covarient_uh

dΩ = Measure(Ω_panel,4*p_fe)
l2(e_u,dΩ)
l2(e_con,dΩ)
l2(e_cov,dΩ)



#### plotting
panel_cfs = [f_panel_cf,grad_covarient_uh,grad_contravarient_uh,e_u,e_con,e_cov,slap_panel_cf]
labels = ["u","sg_cov","sg_con","e_u","e_con","e_cov","slap"]
writevtk_panel(panel_model,panel_cfs,labels)
writevtk_ambient(panel_model,panel_cfs,labels)

### solve the helmholtz problem
e, uh, f_panel_cf = helmholtz_solver(panel_model,f,2)
panel_cfs = [f_panel_cf,uh,f_panel_cf-uh]
labels = ["u","uh","eu"]
writevtk_ambient(panel_model,panel_cfs,labels)
