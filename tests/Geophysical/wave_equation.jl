"""
solve the linearised wave equation in steady form using manufactured solutions
u + ∇ᵧ(φ) = f₁
φ + ∇ᵧ⋅u = f₁
"""

function wave_solver(panel_model,h::Function,vX::Function,p_fe::Int,return_vtk=false)
  lvl = nref(nc(panel_model))
  println("nref = $lvl")

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
  u_proj_cf = panelwise_cellfield(projection_v(vX),Ω_panel,panel_ids)

  sdiv_cf =  panelwise_cellfield(surfdiv(contra_v(vX)),Ω_panel,panel_ids)

  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  sgrad_cf = panelwise_cellfield(sgrad(h),Ω_panel,panel_ids)
  pinvJ_cf = panelwise_cellfield(forward_pinv_jacobian,Ω_panel,panel_ids)

  # manufacture rhs functions
  rhs_scalar = h_cf + sdiv_cf
  rhs_vector = u_proj_cf + sgrad_cf
  rhs_con_vector = pinvJ_cf ⋅ rhs_vector # exact contravariant component


  metric_cf = CellField(analytic_metric,Ω_panel)
  meas_cf = CellField(sqrtg,Ω_panel)
  grad_meas_cf = CellField(grad_meas,Ω_panel)

  biform1((u,p),(v,q)) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
  biform2((u,p),(v,q)) = ∫( (p*q)*meas_cf )dΩ + ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) )  )dΩ

  biformX((u,p),(v,q)) = biform1((u,p),(v,q)) + biform2((u,p),(v,q))
  liformX((v,q)) = ∫( rhs_con_vector⋅(metric_cf⋅v)*meas_cf )dΩ + ∫( (rhs_scalar*q)*meas_cf )dΩ

  op = AffineFEOperator(biformX,liformX,X,Y)
  uh,ph = solve(LUSolver(),op)

  uh_proj = covarient_basis_cf ⋅ uh

  e_u = l2( (u_proj_cf - uh_proj)*meas_cf,dΩ) # error in physical velocity u
  e_p = l2((h_cf - ph)*meas_cf,dΩ) # error in depth

  if return_vtk
    lvl = nref(nc(panel_model))
    cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
    panel_cfs = [h_cf, u_proj_cf, ph, uh_proj, h_cf-ph, u_proj_cf-uh_proj ]
    labels = ["p","u_proj", "ph", "uh_proj", "ep","eu"]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",cellfields=cellfields,append=false,geo_map=cell_geo_map)

  end

  return e_u,e_p
end


function wave_errors(panel_model,h::Function,vX::Function,f::Function,η::Function,p_fe::Int,return_vtk=false)
  e_u,e_p  = wave_solver(panel_model,h,vX,p_fe,return_vtk)
  return e_u,e_p,false
end

function williamson2_convergence_test(solver,n_ref_lvls,return_vtk=false,args...)
  simName = string(solver)[1:end-7]

  println("W2 test")

  for (i,ζ) in enumerate([0.0])
    plot()

    h = panel_to_cartesian(h₀(ζ))
    vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
    f = panel_to_cartesian(f₀(ζ))
    η = panel_to_cartesian(η₀(ζ))

    for p_fe in [1]
      println("p = ", p_fe)
      errs,ns,dxs,slope = convergence_test(solver,n_ref_lvls,h,vX,f,η,p_fe,return_vtk)
      plot_convergence(errs,ns,dxs,slope;
          leginf=["u: p=$p_fe","ϕ: p=$p_fe"],
          colors=[palette(:tab10)[p_fe],palette(:tab10)[p_fe]],
          ls=[:solid, :dot], )


      output = @strdict errs ns dxs slope
      safesave(datadir(dir, ("williamson2_$(simName)_convergence_func_z$(i)_p$p_fe.jld2")), output)
    end
    savefig(plotsdir()*"/williamson2_$(simName)_convergence_func_z$i")
  end

end
