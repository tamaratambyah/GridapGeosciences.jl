#### Linear advection equation
#### solve with SUPG as per Brooks & Hughes 1982 paper
#### Replicate test in Section 5.4 of Rognes2013 paper

vecX(XYZ) = VectorValue(-XYZ[2],XYZ[1],0.0)
u0(XYZ) = exp(-(XYZ[2]^2 + XYZ[3]^2)  )

vX = panel_to_cartesian(tangent_vec(vecX))
u = panel_to_cartesian(u0)


################################################################################
#### Steady with manufactured solutions
################################################################################
function advection_supg_solver(panel_model,u::Function,vX::Function,p_fe::Int,CFL=0.1,return_vtk=false)
  lvl = nref(nc(panel_model))
  println("nref = $lvl")

  panel_ids = get_panel_ids(panel_model)
  degree = 2*(p_fe + 1)

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)

  _rhs(p) = αβ -> u(p)(αβ) + vX(p)(αβ)⋅sgrad(u,p)(αβ)

  v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_cf = panelwise_cellfield(u,Ω_panel,panel_ids)
  rhs_cf = panelwise_cellfield(_rhs,Ω_panel,panel_ids)

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  P = TrialFESpace(Q)

  # hard code RT space as order 1 -- for velocity
  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
  U = TrialFESpace(V)

  vel = interpolate(v_contr_cf,U)
  meas_cf = CellField(sqrtg,Ω_panel)

  # supg stabilisation parameter
  _dx = dx(nc(panel_model))
  dt = _dx*CFL/p_fe
  # τ = (1/dt +  1/_dx )^(-1) # Lima2014
  τ = 0.5*dt

  a_Ω(u,v) = ∫( (u*v)*meas_cf )dΩ + ∫( ((vel⋅∇(u))*v )*meas_cf )dΩ
  a_s(u,v) =  ∫( (u*(vel⋅∇(v)) )*meas_cf )dΩ + ∫( ((vel⋅∇(u))*(vel⋅∇(v)) )*meas_cf )dΩ

  l_Ω(v) = ∫( rhs_cf*v*meas_cf )dΩ
  l_s(v) = ∫( rhs_cf*(vel⋅∇(v))*meas_cf )dΩ

  biform_advection(u,v) = a_Ω(u,v) + τ*a_s(u,v)
  liform_advection(v) = l_Ω(v) + τ*l_s(v)

  op = AffineFEOperator(biform_advection,liform_advection,P,Q)
  uh = solve(LUSolver(),op)

  eu = l2((uh-u_cf)*meas_cf,dΩ)

  if return_vtk
    cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
    labels = ["uh","u","eu"]
    panel_cfs = [uh,u_cf,uh-u_cf]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe", cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  return eu
end


function advection_supg_errors(panel_model,args...)
  e_u  = advection_supg_solver(panel_model,args...)
  return e_u,false,false
end

function advection_supg_convergence_test(n_ref_lvls,u,vX,CFL=0.1,return_vtk=false)
  plot()
  for p_fe in [1,2,3]
    errs,ns,dxs,slope = convergence_test(advection_supg_errors,n_ref_lvls,u,vX,p_fe,CFL,return_vtk)
    plot_convergence(errs,ns,dxs,slope;
        leginf=["u: p=$p_fe"],
        colors=[palette(:tab10)[p_fe]],
        ls=[:solid, :dot], )
  end
  savefig(plotsdir()*"/advection_supg_convergence")

end

n_ref_lvls = 4
CFL = 0.1
advection_supg_convergence_test(n_ref_lvls,u,vX,CFL,true)


################################################################################
#### Transient
################################################################################
function transient_advection_supg(panel_model,u::Function,vX::Function,p_fe::Int,CFL=0.1,return_vtk=false)
  panel_ids = get_panel_ids(panel_model)
  degree = 2*(p_fe + 1)

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)


  v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_cf = panelwise_cellfield(u,Ω_panel,panel_ids)

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  P = TrialFESpace(Q)

  # hard code RT space as order 1 -- for velocity
  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
  U = TrialFESpace(V)

  vel = interpolate(v_contr_cf,U)
  meas_cf = CellField(sqrtg,Ω_panel)

  # supg stabilisation parameter
  _dx = dx(nc(panel_model))
  dt = _dx*CFL/p_fe
  τ = 0.5*dt

  a_mass_Ω(dtu,v) = ∫( (dtu*v)*meas_cf )dΩ
  a_mass_s(dtu,v) = ∫( (dtu*(vel⋅∇(v)))*meas_cf )dΩ
  a_Ω(u,v) = ∫( ((vel⋅∇(u))*v )*meas_cf )dΩ
  a_s(u,v) =  ∫( ((vel⋅∇(u))*(vel⋅∇(v)) )*meas_cf )dΩ

  a_mass(t,dtu,v) = a_mass_Ω(dtu,v) + τ*a_mass_s(dtu,v)
  res(t,u,v) =  a_Ω(u,v) + τ*a_s(u,v)
  jac(t,u,du,v) = a_Ω(du,v) + τ*a_s(du,v)
  jac_t(t,u,dtu,v) = a_mass_Ω(dtu,v) + τ*a_mass_s(dtu,v)
  opT = TransientSemilinearFEOperator(a_mass, res, (jac,jac_t), P, Q, constant_mass=true)

  # solve with SSP RK 3
  uh0 = interpolate(u_cf, P)
  t0, tF = 0.0, 2*π

  solver = RungeKutta(LUSolver(), LUSolver(), dt, :EXRK_SSP_3_3)
  solT = solve(solver, opT, t0, tF, uh0)

  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
  lvl = nref(nc(panel_model))

  println("nlevl = $lvl")
  if return_vtk
    dir = datadir("Transient_advection_nref$lvl")
    !isdir(dir) && mkdir(dir)
    !isdir(plotsdir()) && mkdir(plotsdir())
    writevtk(Ω_panel,dir*"/solT_0.vtu", cellfields=["uh"=>uh0,"v"=>covarient_basis_cf⋅ vel],append=false,geo_map=cell_geo_map)
  end

  ## store errors
  ts = Float64[]
  Es = Float64[]

  push!(ts,0.0)
  push!(Es,0.0)

  for (t,uh) in solT

    println(t)

    eu = l2((uh-uh0)*meas_cf,dΩ)

    push!(ts,t)
    push!(Es,eu)
    if return_vtk
      writevtk(Ω_panel,dir*"/solT_$t.vtu", cellfields=["uh"=>uh],append=false,geo_map=cell_geo_map)
    end
  end

  if return_vtk
    make_pvd(dir,"solT",1)
  end


  plot()
  plot!(ts,Es,lw=3,label="nref = $lvl")
  plot!(xlabel="t",ylabel=L"L2(u_0-u_t)")
  savefig(plotsdir()*"/advection_transient_error_nref$lvl")

  return Es[end]
end

function transient_advection_supg_errors(panel_model,args...)
  e_u  = transient_advection_supg(panel_model,args...)
  return e_u,false,false
end

function transient_advection_supg_convergence_test(n_ref_lvls,u,vX,CFL=0.1,return_vtk=false)

  for p_fe in [1]
    errs,ns,dxs,slope = convergence_test(transient_advection_supg_errors,n_ref_lvls,u,vX,p_fe,CFL,return_vtk)
    plot()
    plot_convergence(errs,ns,dxs,slope;
        leginf=["u: p=$p_fe"],
        colors=[palette(:tab10)[p_fe]],
        ls=[:solid, :dot], )
  end
  savefig(plotsdir()*"/transient_advection_supg_convergence")

end



n_ref_lvls = 4
CFL = 0.2
transient_advection_supg_convergence_test(n_ref_lvls,u,vX,CFL,false)
