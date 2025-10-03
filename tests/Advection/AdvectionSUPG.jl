""" Linear advection equation (material form)
∂ₜu + β⋅ ∇ᵧ(u) = 0
Solve with SUPG as per Brooks & Hughes 1982 paper
Replicate test in Section 5.4 of Rognes2013 paper
"""

using Gridap.Algebra
################################################################################
#### Steady with manufactured solutions
################################################################################
function advection_supg_solver(panel_model,u::Function,vX::Function,p_fe::Int,CFL=0.1,ls=LUSolver(),return_vtk=false)
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
  _dt = _dx*CFL/p_fe
  dt = floor(_dt,sigdigits=1)
  τ = 0.5*dt

  a_Ω(u,v) = ∫( (u*v)*meas_cf )dΩ + ∫( ((vel⋅∇(u))*v )*meas_cf )dΩ
  a_s(u,v) =  ∫( (u*(vel⋅∇(v)) )*meas_cf )dΩ + ∫( ((vel⋅∇(u))*(vel⋅∇(v)) )*meas_cf )dΩ

  l_Ω(v) = ∫( rhs_cf*v*meas_cf )dΩ
  l_s(v) = ∫( rhs_cf*(vel⋅∇(v))*meas_cf )dΩ

  biform_advection(u,v) = a_Ω(u,v) + τ*a_s(u,v)
  liform_advection(v) = l_Ω(v) + τ*l_s(v)

  op = AffineFEOperator(biform_advection,liform_advection,P,Q)

  # uh = solve(ls,op)
  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  uh = FEFunction(P,x)

  eu = l2((uh-u_cf)*meas_cf,dΩ)

  if return_vtk
    cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
    labels = ["uh","u","eu"]
    panel_cfs = [uh,u_cf,uh-u_cf]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe", cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  return eu,false,false
end


function advection_supg_convergence_test(dir,u::Function,vX::Function,n_ref_lvls=4,ps=[1],CFL=0.1,ls=LUSolver(),return_vtk=false)
  println("serial advection supg test")

  models  = get_refined_models(n_ref_lvls)

  simName = "advection_supg_convergence_func"

  errors = Vector{Vector{Float64}}(undef,length(ps))
  ns = Vector{Vector{Float64}}(undef,length(ps))
  dxs = Vector{Vector{Float64}}(undef,length(ps))
  slopes = Vector{Float64}(undef,length(ps))

  for (i,p_fe) in enumerate(ps)
    println("p_fe = $p_fe")
    errors[i],ns[i],dxs[i],slopes[i] = h_convergence_test(models,advection_supg_solver,u,vX,p_fe,CFL,ls,return_vtk)
  end

  print_convergence_results(errors,ns,dxs,slopes,ps)
  output = @strdict errors ns dxs slopes ps

  safesave(datadir(dir, ("$simName.jld2")), output)

  plot_convergence_from_saved(dir,simName)


end



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
  _dt = _dx*CFL/p_fe
  dt = floor(_dt,sigdigits=1)
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

  return Es[end],false,false
end

# function transient_advection_supg_errors(panel_model,args...)
#   e_u  = transient_advection_supg(panel_model,args...)
#   return e_u,false,false
# end

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
