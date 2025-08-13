## make models
cube_model = coarse_cube_model(π/4,6)
cube_model = Gridap.Adaptivity.refine(cube_model)
cube_model = Gridap.Adaptivity.refine(cube_model)

panel_model,panel_ids = parametric_model(cube_model)

sphere_model = ambient_model(panel_model,panel_ids)

Ω_panel = Triangulation(panel_model)
Ω_sphere = Triangulation(sphere_model)

l2(e,dΩ) = sum(∫( e⋅e )dΩ)


################################################################################
#### On panel
################################################################################
### Panel cell fields
f_panel_cf = panelwise_cellfield(f_sin,Ω_panel,panel_ids)
covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
contravarient_basis_cf = panelwise_cellfield(contravariant_basis,Ω_panel,panel_ids)
inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)

### analytic gradient -- need to define new function to trigger automatric differentiation
gradf(p) = αβ -> gradient(f_sin(p))(αβ)
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

e_con = grad_contravarient-grad_contravarient_uh
e_cov = grad_covarient-grad_covarient_uh

dΩ = Measure(Ω_panel,4*p_fe)
l2(e_con,dΩ)
l2(e_cov,dΩ)


#### plotting
panel_cfs = [f_uh,grad_covarient,grad_contravarient,e_con,e_cov]
labels = ["sinu","sg_cov","sg_con","e_con","e_cov"]
writevtk_panel(panel_model,panel_cfs,labels,panel_ids)
writevtk_ambient(ambient_model,panel_cfs,labels,panel_ids)
