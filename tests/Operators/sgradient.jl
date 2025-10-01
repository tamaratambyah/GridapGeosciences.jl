"""
test the computation of surface gradient
"""

function sgrad(panel_model,func::Function,p_fe::Int)
  lvl = nref(nc(panel_model))
  println("nref = $lvl")

  Ω_panel = Triangulation(panel_model)
  panel_ids = get_panel_ids(panel_model)

  f_panel_cf = panelwise_cellfield(func,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  contravarient_basis_cf = panelwise_cellfield(contravariant_basis,Ω_panel,panel_ids)
  inv_metric_cf = CellField(analytic_inv_metric,Ω_panel)

  ### analytic gradient -- need to define new function to trigger automatric differentiation
  gradf(p) = αβ -> gradient(func(p))(αβ)
  gradf_cf = panelwise_cellfield(gradf,Ω_panel,panel_ids)

  grad_covarient =  contravarient_basis_cf ⋅ gradf_cf
  grad_contravarient = covarient_basis_cf ⋅  (inv_metric_cf ⋅ gradf_cf)

  ### interpolate f into FE space, and then differentiate
  V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  U = TrialFESpace(V)

  f_uh = interpolate(f_panel_cf,U)

  grad_covarient_uh =  contravarient_basis_cf ⋅ gradient(f_uh)
  grad_contravarient_uh = covarient_basis_cf ⋅  (inv_metric_cf ⋅ gradient(f_uh))

  return grad_covarient, grad_contravarient, grad_covarient_uh,grad_contravarient_uh

end


function sgrad_errors(panel_model,func::Function,p_fe::Int)

  grad_covarient, grad_contravarient, grad_covarient_uh,grad_contravarient_uh = sgrad(panel_model,func,p_fe)

  e_con = grad_contravarient-grad_contravarient_uh
  e_cov = grad_covarient-grad_covarient_uh

  dΩ = Measure(Triangulation(panel_model),4*p_fe)
  e1 = l2(e_con,dΩ)
  e2 = l2(e_cov,dΩ)

  return e1,e2,false
end

function sgrad_convergence_test(analytic_funcs,n_ref_lvls)

  for (key, val) in analytic_funcs
    plot()
    for p_fe in [1,2,3]
      errs,ns,dxs,slope = convergence_test(sgrad_errors,n_ref_lvls,val,p_fe)
      plot_convergence(errs,ns,dxs,slope;leginf=["p=$(p_fe)_con","p=$(p_fe)_cov"],colors=palette(:tab10))
    end
    savefig(plotsdir()*"/sgrad_convergence_func_$(key)")
  end
end
