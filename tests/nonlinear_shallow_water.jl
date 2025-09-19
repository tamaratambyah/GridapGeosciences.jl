function nonlinear_shallow_water_solver(panel_model,h::Function,vX::Function,f::Function,η::Function,
    p_fe::Int,return_vtk=false,check_geo_balance=false)

  lvl = nref(nc(panel_model))
  println("Refinement level: $lvl")

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  R = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1)
  H = TrialFESpace(R)

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)


  ## initial conditions
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_contra_h = interpolate(u_contra_cf,U)
  u_proj_h = covarient_basis_cf ⋅ u_contra_h

  h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
  h_h = interpolate(h_cf,P)

  cor_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  gravity = 1.0

  # absolute vorticity
  η_cf = panelwise_cellfield(η,Ω_panel,panel_ids)
  η_h = interpolate(η_cf,H)

  # mectrics required in weak forms
  detg_cf = CellField(detg,Ω_panel)
  metric_cf = CellField(analytic_metric,Ω_panel)
  inv_metric_cf = CellField(analytic_inv_metric,Ω_panel)
  meas_cf = CellField(sqrtg,Ω_panel)
  grad_meas_cf = CellField(grad_meas,Ω_panel)


  #### DIAGNOSTIC VARIABLES
  # mass flux
  biformF(F,v) = ∫( (F⋅ (metric_cf⋅v))*meas_cf )dΩ
  liformF(v) = ∫( h_h*(u_contra_h⋅(metric_cf⋅v))*meas_cf   )dΩ
  op = AffineFEOperator(biformF,liformF,U,V)
  Fh = solve(LUSolver(),op)

  # Bernoulli potential
  biformΦ(Φ,r) = ∫( Φ*r*meas_cf  )dΩ
  liformΦ(r) = ∫( gravity*h_h*r*meas_cf  )dΩ + ∫( 0.5*( u_contra_h ⋅(metric_cf⋅u_contra_h) )r*meas_cf  )dΩ
  op = AffineFEOperator(biformΦ,liformΦ,P,Q)
  Φh = solve(LUSolver(),op)

  # vorticity
  perp_matrix_cf = CellField(analytic_perp_matrix,Ω_panel)
  biformq(q,r) = ∫( q*h_h*r*meas_cf  )dΩ
  liformq(r) = ∫( cor_cf*r*meas_cf  )dΩ + ∫( (perp_matrix_cf⋅u_contra_h)⋅∇(r)  )dΩ
  op = AffineFEOperator(biformq,liformq,H,R)
  qh = solve(LUSolver(),op)

  e_η = l2((η_h - qh*h_h )*meas_cf,dΩ)

  #### PROGNOSTIC VARIABLES

  # equation for depth:
  rhs_h = h_h + 1/meas_cf*( Fh⋅grad_meas_cf + meas_cf*(∇⋅Fh)   )
  biform_p(p,r) = ∫( (p*r)*meas_cf )dΩ
  liform_p(r) = ( ∫( (rhs_h*r)*meas_cf )dΩ
                - ∫( r*(Fh⋅grad_meas_cf + meas_cf*(∇⋅Fh) )  )dΩ
                )
  op = AffineFEOperator(biform_p,liform_p,P,Q)
  ph = solve(LUSolver(),op)
  e_p = l2((h_cf - ph)*meas_cf,dΩ) # error in depth


  # equation for velocity
  Fperph = 1/meas_cf*(perp_matrix_cf⋅Fh)
  rhs_u = u_contra_h + qh*Fperph + (  inv_metric_cf⋅gradient(Φh) )

  # check geostropohic balance
  geo_balance = qh*Fperph + (  inv_metric_cf⋅gradient(Φh) )
  e_geo_balance = sum(∫( geo_balance )dΩ)
  if check_geo_balance
    @check vector_length(e_geo_balance) < 1e-12
    println("Global geostropohic balance error: $e_geo_balance")
  end

  # solve for velocity
  biform_u(u,v) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ
  liform_u(v) = ( ∫( rhs_u⋅(metric_cf⋅v)*meas_cf )dΩ
                  + ∫( Φh*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
                  - ∫( qh*( (perp_matrix_cf⋅Fh) ⋅(metric_cf ⋅v))   )dΩ
                    )
  op = AffineFEOperator(biform_u,liform_u,U,V)
  uh = solve(LUSolver(),op)

  uh_proj = covarient_basis_cf ⋅ uh

  # e_u = l2( ( uh-u_contra_h  )*meas_cf,dΩ  )
  e_u = l2( (u_proj_h - uh_proj)*meas_cf,dΩ) # error in physical velocity u

  if return_vtk
    lvl = nref(nc(panel_model))
    cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
    panel_cfs = [ph, h_cf,ph-h_cf,
                uh_proj, u_proj_h, uh_proj-u_proj_h,
                 ]
    labels = ["ph","p","ep",
              "uh","u","eu",
                ]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$lvl",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  println(e_u, "; ", e_p, "; ",  e_η)

  return e_u,  e_p, e_η, e_geo_balance

end



function nonlinear_shallow_water_errors(panel_model,h::Function,vX::Function,f::Function,η::Function,p_fe::Int,
    return_vtk=false,check_geo_balance=false)
  e_u,e_p,e_η,e_geo_balance  = nonlinear_shallow_water_solver(panel_model,h,vX,f,η,p_fe,return_vtk,check_geo_balance)
  return e_u,e_p,e_η
end


function williamson2_convergence_test(n_ref_lvls,return_vtk=false,args...)

  for (i,ζ) in enumerate([0, π/2 ])

    println("ζ = $ζ")

    plot()

    h = panel_to_cartesian(h₀(ζ))
    vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
    f = panel_to_cartesian(f₀(ζ))
    η = panel_to_cartesian(η₀(ζ))


    for p_fe in [1]
      errs,ns,dxs,slope = convergence_test(nonlinear_shallow_water_errors,n_ref_lvls,h,vX,f,η,p_fe,return_vtk,true)
      plot_convergence(errs,ns,dxs,slope;
          leginf=["u: p=$p_fe","ϕ: p=$p_fe", "η: p=$p_fe"],
          colors=[palette(:tab10)[p_fe],palette(:tab10)[p_fe],palette(:tab10)[p_fe]],
          ls=[:solid, :dot,:dashdot], )
    end
    savefig(plotsdir()*"/williamson2_NL_sw_convergence_func_z$i")
  end

end


## Williamson2 convergence test
n_ref_lvls = 4
williamson2_convergence_test(n_ref_lvls,true)
