### convergence tests
function interpolation(panel_model,func::Function,p_fe::Int)

  Ω_panel = Triangulation(panel_model)
  panel_ids = get_panel_ids(panel_model)

  f_panel_cf = panelwise_cellfield(func,Ω_panel,panel_ids)

  ### interpolate f into FE space, and then differentiate
  V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  U = TrialFESpace(V)

  f_uh = interpolate(f_panel_cf,U)

  return f_panel_cf, f_uh

end


function interpolation_errors(panel_model,func::Function,p_fe::Int)

  f_panel_cf, f_uh = interpolation(panel_model,func,p_fe)

  e  = f_panel_cf - f_uh
  dΩ = Measure(Triangulation(panel_model),4*p_fe)
  e1 = l2(e,dΩ)
  e2 = false

  return e1,e2,false
end

function interpolation_convergence_test(analytic_funcs,n_ref_lvls)

  for (key, val) in analytic_funcs
    plot()
    for p_fe in [1,2,3]
      errs,ns,dxs,slope = convergence_test(interpolation_errors,n_ref_lvls,val,p_fe)
      plot_convergence(errs,ns,dxs,slope;leginf=["p=$p_fe"],colors=[palette(:tab10)[p_fe]])
    end
    savefig(plotsdir()*"/interpolation_convergence_func_$(key)")
  end

end


analytic_funcs = Dict{Symbol,Any}()
analytic_funcs[:sin] = f_sin
analytic_funcs[:XYZ] = f_XYZ

n_ref_lvls = 4

interpolation_convergence_test(analytic_funcs,n_ref_lvls)
