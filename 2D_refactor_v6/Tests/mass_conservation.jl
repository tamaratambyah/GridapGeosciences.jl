
################################################################################
############# DARCY
################################################################################
f =  f_XYZ
f_panel_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  sigma_cf = panelwise_cellfield(contr_gradf(f),Ω_panel,panel_ids)
  sdiv_cf =  panelwise_cellfield(surfdiv(contr_gradf(f)),Ω_panel,panel_ids)
  slap_panel_cf =  panelwise_cellfield(surflap(f),Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  gradu_cf = covarient_basis_cf ⋅ sigma_cf

  rhs_cf = -slap_panel_cf

  V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  U = TrialFESpace(V)

  T = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  S = TrialFESpace(T)

  Y = MultiFieldFESpace([V, T])
  X = MultiFieldFESpace([U, S])

  metric_cf = CellField(analytic_metric,Ω_panel)
  meas_cf = CellField(sqrtg,Ω_panel)
  grad_meas_cf = CellField(grad_meas,Ω_panel)

  biform1((u,s),(v,t)) = ∫( (s⋅ (metric_cf⋅t))*meas_cf )dΩ + ∫( u*(t⋅grad_meas_cf + meas_cf*(∇⋅t) ) )dΩ
  biform2((u,s),(v,t)) = ∫( -( v*(s⋅grad_meas_cf + meas_cf*(∇⋅s) ) ))dΩ

  biformX((u,s),(v,t)) = biform1((u,s),(v,t)) + biform2((u,s),(v,t))
  liformX((v,t)) = ∫( (rhs_cf*v)*meas_cf )dΩ

  op = AffineFEOperator(biformX,liformX,X,Y)
  uh,sh = solve(LUSolver(),op)
  graduh = covarient_basis_cf ⋅sh

  e_u = l2( (f_panel_cf - uh)*meas_cf,dΩ) # error in scalar u
  e_s = l2((sigma_cf - sh)*meas_cf,dΩ) # error in contra compons of grad u
  e_gradu = l2((gradu_cf - graduh)*meas_cf,dΩ) # error in grad u = physical sigma


sum(∫( 1/meas_cf *  divergence(meas_cf*sh) )dΩ)
# sum(∫(  1/meas_cf * ( uh⋅grad_meas_cf + meas_cf*(divergence(uh))) )dΩ)
sum(∫( divergence(sh) )dΩ)
