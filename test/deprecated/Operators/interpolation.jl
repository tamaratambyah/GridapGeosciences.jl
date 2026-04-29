"""
test the interpolation of function into FE space
"""

### convergence tests
function interpolation(panel_model,p_fe::Int,dir::String,func::Function)
  lvl = nref(nc(panel_model))
  println("p_fe = $(p_fe); nref = $lvl")

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*p_fe+1)
  panel_ids = get_panel_ids(panel_model)

  f_panel_cf = ParametricCellField(func,Ω_panel,panel_ids)

  ### interpolate f into FE space, and then differentiate
  V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  U = TrialFESpace(V)

  f_uh = interpolate(f_panel_cf,U)
  e  = l2(f_panel_cf - f_uh,dΩ)
  return e,false,false

end



function interpolation_convergence_test(ranks::AbstractArray,nprocs::Int,analytic_funcs,n_ref_lvls=4,ps=[1,2,3])
  # serial models
  models  = get_refined_models(n_ref_lvls)
  dir = datadir("InterpolationConvergence")
  !isdir(dir) && mkdir(dir)

  for (key, val) in analytic_funcs
    _dir = dir*"/func_$(key)"
    !isdir(_dir) && mkdir(_dir)
    p_convergence_test(ranks,ps,models,interpolation,_dir,val)
    plot_convergence_from_saved(_dir,"convergence",["p"])
  end

end
