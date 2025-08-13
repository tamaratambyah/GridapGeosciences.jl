function helmholtz_solver(panel_model::ParametricDiscreteModel,f::Function,p_fe::Int)
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*p_fe+1)

  V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  U = TrialFESpace(V)

  f_panel_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  inv_metric_cf = CellField(analytic_inv_metric,Ω_panel)
  meas_cf = CellField(sqrtg,Ω_panel)
  slap_panel_cf =  panelwise_cellfield(surflap(f),Ω_panel,panel_ids)

  rhs_cf = f_panel_cf + slap_panel_cf

  poisson_biform(u,v) = ∫(u*v*meas_cf)dΩ -  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_cf )dΩ
  poisson_liform(v) = ∫(  (rhs_cf*v)*meas_cf )dΩ
  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
  uh = solve(LUSolver(),op)

  e = l2(f_panel_cf-uh,dΩ)
  return e, uh, f_panel_cf
end


function helmholtz_errors(panel_model::ParametricDiscreteModel,func::Function,p_fe::Int)
  e,  = helmholtz_solver(panel_model,func,p_fe)
  return e,false
end

function helmholtz_convergence_test(analytic_funcs,n_ref_lvls)

  for (key, val) in analytic_funcs
    plot()
    for p_fe in [1,2,3]
      errs,ns,dxs,slope = convergence_test(helmholtz_errors,n_ref_lvls,val,p_fe)
      plot_convergence(errs,ns,dxs,slope;leginf=["p=$p_fe"],colors=[palette(:tab10)[p_fe]])
    end
    savefig(plotsdir()*"/helmholtz_convergence_func_$(key)")
  end

end


analytic_funcs = Dict{Symbol,Any}()
analytic_funcs[:sin] = f_sin
analytic_funcs[:XYZ] = f_XYZ

n_ref_lvls = 4

helmholtz_convergence_test(analytic_funcs,n_ref_lvls)

#### mapped functions
mapped_funcs = Dict{Symbol,Any}()
mapped_funcs[:sin] = panel_to_latlon(fθϕ)
mapped_funcs[:XYZ] = panel_to_cartesian(fX)

n_ref_lvls = 4
helmholtz_convergence_test(mapped_funcs,n_ref_lvls)



#### Williamson 2
williamson_funcs = Dict{Symbol,Any}()
williamson_funcs[:z1] = panel_to_latlon(fWilliamson(0))
williamson_funcs[:z2] = panel_to_latlon(fWilliamson(0.05))
williamson_funcs[:z3] = panel_to_cartesian(fWilliamson(π/2-0.05))
williamson_funcs[:z4] = panel_to_cartesian(fWilliamson(π/2))

n_ref_lvls = 4
helmholtz_convergence_test(williamson_funcs,n_ref_lvls)
