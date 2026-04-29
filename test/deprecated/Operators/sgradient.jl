"""
test the computation of surface gradient
"""

function sgrad(panel_model,p_fe::Int,dir::String,func::Function)
  lvl = nref(nc(panel_model))
  println("p_fe = $(p_fe); nref = $lvl")

  Ω_panel = Triangulation(panel_model)
  panel_ids = get_panel_ids(panel_model)

  f_panel_cf = ParametricCellField(func,Ω_panel,panel_ids)
  covariant_basis_cf = ParametricCellField(covariant_basis,Ω_panel,panel_ids)
  contravarient_basis_cf = ParametricCellField(contravariant_basis,Ω_panel,panel_ids)
  inv_metric_cf = ParametricCellField(inv_metric,Ω_panel,panel_ids)

  ### analytic gradient -- need to define new function to trigger automatric differentiation
  gradf(p) = αβ -> gradient(func(p))(αβ)
  gradf_cf = ParametricCellField(gradf,Ω_panel,panel_ids)

  grad_covarient =  contravarient_basis_cf ⋅ gradf_cf
  grad_contravarient = covariant_basis_cf ⋅  (inv_metric_cf ⋅ gradf_cf)

  ### interpolate f into FE space, and then differentiate
  V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  U = TrialFESpace(V)

  f_uh = interpolate(f_panel_cf,U)

  grad_covarient_uh =  contravarient_basis_cf ⋅ gradient(f_uh)
  grad_contravarient_uh = covariant_basis_cf ⋅  (inv_metric_cf ⋅ gradient(f_uh))

  e_con = grad_contravarient-grad_contravarient_uh
  e_cov = grad_covarient-grad_covarient_uh

  dΩ = Measure(Ω_panel,4*p_fe)
  e1 = l2(e_con,dΩ)
  e2 = l2(e_cov,dΩ)

  return e1,e2,false

end


function sgrad_convergence_test(ranks::AbstractArray,nprocs::Int,analytic_funcs,n_ref_lvls=4,ps=[1,2,3])
  # serial models
  models  = get_refined_models(n_ref_lvls)
  dir = datadir("SgradConvergence")
  !isdir(dir) && mkdir(dir)

  for (key, val) in analytic_funcs
    _dir = dir*"/func_$(key)"
    !isdir(_dir) && mkdir(_dir)
    p_convergence_test(ranks,ps,models,sgrad,_dir,val)
    plot_convergence_from_saved(_dir,"convergence",["con","cov"])
  end

end
