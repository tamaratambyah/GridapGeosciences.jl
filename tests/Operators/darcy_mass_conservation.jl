"""
test mass conservation by solving
u + ∇ᵧ(φ) = u₀
∇ᵧ⋅u = 0
for arbitray vector field u₀
"""


function mass_conservation(panel_model,func::Function,p_fe::Int,scalar_field::Bool,return_vtk=false)
  lvl = nref(nc(panel_model))
  println("nref = $lvl")

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  vec_contra_cf = if scalar_field
    panelwise_cellfield(contr_gradf(func),Ω_panel,panel_ids)
  else
    panelwise_cellfield(contra_v(func),Ω_panel,panel_ids)
  end

  metric_cf = CellField(analytic_metric,Ω_panel)
  meas_cf = CellField(sqrtg,Ω_panel)
  grad_meas_cf = CellField(grad_meas,Ω_panel)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  biform1((u,p),(v,q)) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
  biform2((u,p),(v,q)) =  ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) )  )dΩ

  biformX((u,p),(v,q)) = biform1((u,p),(v,q)) + biform2((u,p),(v,q))
  liformX((v,q)) = ∫( vec_contra_cf⋅(metric_cf⋅v)*meas_cf )dΩ

  op = AffineFEOperator(biformX,liformX,X,Y)
  uh,ph = solve(LUSolver(),op)

  # mass conservation errors
  s_div = sum(∫(  divergence(meas_cf*uh) )dΩ)
  s_div0 = sum(∫(  divergence(meas_cf*vec_contra_cf) )dΩ)
  panel_div = sum(∫(  divergence(uh) )dΩ)

  if return_vtk
    lvl = nref(nc(panel_model))
    cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)

    u_proj = covarient_basis_cf ⋅ vec_contra_cf
    u_projh = covarient_basis_cf ⋅ uh

    panel_cfs = [ u_proj, u_projh, ph]
    labels = ["u0", "u_projh", "p"]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  return s_div, s_div0, panel_div

end


function mass_conservation_errors(panel_model,func::Function,p_fe::Int,scalar_field::Bool,return_vtk=false)
  s_div, s_div0, = mass_conservation(panel_model,func,p_fe,scalar_field,return_vtk )
  println("initial divergence: $s_div0")
  return abs(s_div),false,false
end

function mass_conservation_convergence_test(analytic_funcs,n_ref_lvls,scalar_field::Bool,return_vtk=false)

  for (key, val) in analytic_funcs
    plot()
    for p_fe in [1,2]
      errs,ns,dxs,slope = convergence_test(mass_conservation_errors,n_ref_lvls,val,p_fe,scalar_field,return_vtk)
      plot_error(ns,errs;
          leginf=["dM: p=$p_fe"],
          colors=[palette(:tab10)[p_fe]],
          ls=[:solid, :dot], )
      plot!(yscale=:log10,framestyle=:box,
          xscale=:log10,xlabel="n cells",ylabel="dM"
          )
    end
    savefig(plotsdir()*"/darcy_mass_conservation_convergence_func_$(key)")
  end

end
