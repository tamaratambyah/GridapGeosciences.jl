
function mixed_helmholtz_solver(panel_model,f::Function,p_fe::Int,return_vtk=false)

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  f_panel_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  sigma_cf = panelwise_cellfield(contr_gradf(f),Ω_panel,panel_ids)
  sdiv_cf =  panelwise_cellfield(surfdiv(contr_gradf(f)),Ω_panel,panel_ids)
  slap_panel_cf =  panelwise_cellfield(surflap(f),Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  gradu_cf = covarient_basis_cf ⋅ sigma_cf

  rhs_cf = f_panel_cf + sdiv_cf

  println("Check sdiv(sgrad) against slap: ", l2(sdiv_cf-slap_panel_cf,dΩ) )


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
  biform2((u,s),(v,t)) = ∫( (u*v)*meas_cf )dΩ + ∫( v*(s⋅grad_meas_cf + meas_cf*(∇⋅s) ) )dΩ

  biformX((u,s),(v,t)) = biform1((u,s),(v,t)) + biform2((u,s),(v,t))
  liformX((v,t)) = ∫( (rhs_cf*v)*meas_cf )dΩ

  op = AffineFEOperator(biformX,liformX,X,Y)
  uh,sh = solve(LUSolver(),op)
  graduh = covarient_basis_cf ⋅sh


  e_u = l2( (f_panel_cf - uh)*meas_cf,dΩ) # error in scalar u
  e_s = l2((sigma_cf - sh)*meas_cf,dΩ) # error in contra compons of grad u
  e_gradu = l2((gradu_cf - graduh)*meas_cf,dΩ) # error in grad u = physical sigma

  if return_vtk
    lvl = nref(nc(panel_model))
    cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)

    panel_cfs = [uh, sh, graduh, f_panel_cf, gradu_cf, rhs_cf, f_panel_cf - uh, sigma_cf - sh, gradu_cf - graduh  ]
    labels = ["uh","sh","graduh", "u_ex", "gradu","rhs", "eu", "es", "egradu"]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$lvl",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  return e_u,e_s,e_gradu

end



function mixed_helmholtz_errors(panel_model,func::Function,p_fe::Int,return_vtk=false)
  e_u,e_s,e_gradu  = mixed_helmholtz_solver(panel_model,func,p_fe,return_vtk)
  return e_u,e_gradu,false
end

function mixed_helmholtz_convergence_test(analytic_funcs,n_ref_lvls,return_vtk=false)

  for (key, val) in analytic_funcs
    plot()
    for p_fe in [1,2]
      errs,ns,dxs,slope = convergence_test(mixed_helmholtz_errors,n_ref_lvls,val,p_fe,return_vtk)
      plot_convergence(errs,ns,dxs,slope;
          leginf=["u: p=$p_fe","s: p=$p_fe"],
          colors=[palette(:tab10)[p_fe],palette(:tab10)[p_fe]],
          ls=[:solid, :dot], )
    end
    savefig(plotsdir()*"/mixed_helmholtz_convergence_func_$(key)")
  end

end


analytic_funcs = Dict{Symbol,Any}()
analytic_funcs[:sin] = f_sin
analytic_funcs[:XYZ] = f_XYZ

n_ref_lvls = 4

mixed_helmholtz_convergence_test(analytic_funcs,n_ref_lvls,false)

#### mapped functions
mapped_funcs = Dict{Symbol,Any}()
mapped_funcs[:sin] = panel_to_latlon(fθϕ)
mapped_funcs[:XYZ] = panel_to_cartesian(fX)

n_ref_lvls = 4
mixed_helmholtz_convergence_test(mapped_funcs,n_ref_lvls)



#### Williamson 2
williamson_funcs = Dict{Symbol,Any}()
williamson_funcs[:z1] = panel_to_latlon(fWilliamson(0))
williamson_funcs[:z2] = panel_to_latlon(fWilliamson(0.05))
williamson_funcs[:z3] = panel_to_cartesian(fWilliamson(π/2-0.05))
williamson_funcs[:z4] = panel_to_cartesian(fWilliamson(π/2))

n_ref_lvls = 4
mixed_helmholtz_convergence_test(williamson_funcs,n_ref_lvls)
