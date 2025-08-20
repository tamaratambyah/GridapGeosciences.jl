function wave_solver(panel_model,h::Function,vX::Function,p_fe::Int,return_vtk=false)

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
    writevtk(Ω_panel,dir*"/ambient_model_nref$lvl",cellfields=cellfields,append=false,geo_map=cell_geo_map)

  end

  return e_u,e_p
end


function wave_errors(panel_model,h::Function,vX::Function,p_fe::Int,return_vtk=false)
  e_u,e_p  = wave_solver(panel_model,h,vX,p_fe,return_vtk)
  return e_u,e_p
end

function williamson2_convergence_test(n_ref_lvls,return_vtk=false,args...)

  for (i,ζ) in enumerate([π/2])
    plot()

    h = panel_to_latlon(hWilliamson(ζ,args...))
    vecX = vec_cartesian_to_latlon(vWilliamson(ζ,args...))
    vX = panel_to_cartesian(tangent_vec(vecX))

    for p_fe in [1]
      errs,ns,dxs,slope = convergence_test(wave_errors,n_ref_lvls,h,vX,p_fe,return_vtk)
      plot_convergence(errs,ns,dxs,slope;
          leginf=["u: p=$p_fe","ϕ: p=$p_fe"],
          colors=[palette(:tab10)[p_fe],palette(:tab10)[p_fe]],
          ls=[:solid, :dot], )
    end
    savefig(plotsdir()*"/williamson2_wave_convergence_func_z$i")
  end

end


vWilliamson(ζ,u0,ω) = θϕ -> - u0*VectorValue( cos(θϕ[2])*cos(ζ) + cos(θϕ[1])*sin(θϕ[2])*sin(ζ),
                                      -sin(θϕ[1])*sin(ζ) )

hWilliamson(ζ,u0,ω) = θϕ -> 1 - (ω*u0 + 0.5*u0^2)*( -cos(θϕ[1])*cos(θϕ[2])*sin(ζ) +  sin(θϕ[2])*cos(ζ) )^2

ω = 1e-5
u0 = 0.1
n_ref_lvls = 4
williamson2_convergence_test(n_ref_lvls,true,u0,ω)
