""" Linear advection equation (flux form)
∂ₜu + ∇ᵧ⋅(βu) = 0
Solve with dG upwinding as per Brezzi 2004 paper
Replicate test in Section 5.4 of Rognes2013 paper
"""

using Gridap.Geometry, Gridap.Algebra

function my_mean( Bu_n::SkeletonPair)
  plus  = ( Bu_n.plus)
  minus = ( Bu_n.minus)
  0.5*( plus - minus  )
end

function _my_mean(j::SkeletonPair,vel::CellField,u::CellField)
  0.5*( (j.plus⋅vel.plus)*u.plus + (j.minus⋅vel.minus)*u.minus )
end

function my_jump(j::SkeletonPair,ginv::SkeletonPair,n::SkeletonPair,u::CellField)
  u.plus*(j.plus⋅(ginv.plus⋅n.plus) ) + u.minus*(j.minus⋅(ginv.minus⋅n.minus) )
end


################################################################################
#### Steady with manufactured solutions
################################################################################
function advection_dg_solver(panel_model,u::Function,vX::Function,uvX::Function,p_fe::Int,ls=LUSolver(),return_vtk=false)
  lvl = nref(nc(panel_model))
  println("nref = $lvl")

  panel_ids = get_panel_ids(panel_model)
  degree = 2*(p_fe + 1)

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)

  Λ = SkeletonTriangulation(panel_model)
  dΛ = Measure(Λ,degree)
  n_Λ = get_normal_vector(Λ)

  _rhs(p) = αβ -> u(p)(αβ) + surfdiv(contra_v(uvX))(p)(αβ)

  v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_cf = panelwise_cellfield(u,Ω_panel,panel_ids)
  rhs_cf = panelwise_cellfield(_rhs,Ω_panel,panel_ids)

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  # hard code RT space as order 1 -- for velocity
  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
  U = TrialFESpace(V)

  vel = interpolate(v_contr_cf,U)
  meas_cf = CellField(sqrtg,Ω_panel)

  a_Ω(u,v) = ∫( (u*v)*meas_cf )dΩ - ∫( (u*(∇(v)⋅vel) )*meas_cf )dΩ

  ### volume stabilisation term
  # a_s1(u,v) = ∫( my_mean((vel*u)⋅n_Λ)*jump(v)*meas_cf   )dΛ

  jac_cf = panelwise_cellfield(forward_jacobian,Λ)
  ginv_cf = panelwise_cellfield(_analytic_inv_metric,Λ)
  a_s1(u,v) = ∫( _my_mean(jac_cf,vel,u)⋅my_jump(jac_cf,ginv_cf,n_Λ,v)*meas_cf   )dΛ


  ### upwinding stabilisation term
  upwind = abs((vel⋅ n_Λ).plus)
  # a_s2(u,v) = ∫(  0.5*(upwind)*jump(u)*jump(v)*meas_cf   )dΛ

  cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
  cell_normal = get_facet_normal(Λ,cell_geo_map)
  n = get_normal_vector(Λ,cell_normal)
  a_s2(u,v) = ∫(  0.5*(upwind)*jump(u*n)⋅jump(v*n)*meas_cf   )dΛ


  biform_advection(p,q) =  a_Ω(p,q) + a_s1(p,q) + a_s2(p,q)
  liform_advection(q) = ∫( (rhs_cf*q)*meas_cf )dΩ

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
    lvl = nref(nc(panel_model))
    cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
    labels = ["uh","u","eu"]
    panel_cfs = [uh,u_cf,uh-u_cf]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe", cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end
  return eu,false,false
end



function advection_dg_convergence_test(dir,u::Function,vX::Function,uvX::Function,n_ref_lvls=4,ps=[1],ls=LUSolver(),return_vtk=false)
  println("serial advection DG test")

  models  = get_refined_models(n_ref_lvls)

  simName = "advection_dg_convergence_func"

  errors = Vector{Vector{Float64}}(undef,length(ps))
  ns = Vector{Vector{Float64}}(undef,length(ps))
  dxs = Vector{Vector{Float64}}(undef,length(ps))
  slopes = Vector{Float64}(undef,length(ps))

  for (i,p_fe) in enumerate(ps)
    println("p_fe = $p_fe")
    errors[i],ns[i],dxs[i],slopes[i] = h_convergence_test(models,advection_dg_solver,u,vX,uvX,p_fe,ls,return_vtk)
  end

  print_convergence_results(errors,ns,dxs,slopes,ps)
  output = @strdict errors ns dxs slopes ps

  safesave(datadir(dir, ("$simName.jld2")), output)

  plot_convergence_from_saved(dir,simName)


end


################################################################################
#### Transient
################################################################################
function transient_advection_dg(panel_model,u::Function,vX::Function,p_fe::Int,CFL=0.1,return_vtk=false)
  lvl = nref(nc(panel_model))
  println("nref = $lvl")

  panel_ids = get_panel_ids(panel_model)
  degree = 2*(p_fe + 1)

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)

  Λ = SkeletonTriangulation(panel_model)
  dΛ = Measure(Λ,degree)
  n_Λ = get_normal_vector(Λ)

  v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_cf = panelwise_cellfield(u,Ω_panel,panel_ids)

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TransientTrialFESpace(Q)

  # hard code RT space as order 1 -- for velocity
  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
  U = TrialFESpace(V)

  # initial conditions
  vel = interpolate(v_contr_cf,U)
  uh0 = interpolate(u_cf, P(0.0))

  meas_cf = CellField(sqrtg,Ω_panel)

  ## weak form
  a_mass(t,dtu,v) = ∫( (dtu*v)*meas_cf )dΩ

  a_Ω(u,v) =   ∫( -(u*(∇(v)⋅vel) )*meas_cf )dΩ
  a_s1(u,v) = ∫( my_mean((vel*u)⋅n_Λ)*jump(v)*meas_cf   )dΛ

  upwind = abs((vel⋅ n_Λ).plus)/2
  a_s2(u,v) = ∫(  upwind*jump(u)*jump(v)*meas_cf   )dΛ

  res(t,u,v) =  a_Ω(u,v) + a_s1(u,v) + a_s2(u,v)
  jac(t,u,du,v) = a_Ω(du,v) + a_s1(du,v) + a_s2(du,v)
  jac_t(t,u,dtu,v) = ∫( (dtu*v)*meas_cf )dΩ
  opT = TransientSemilinearFEOperator(a_mass, res, (jac,jac_t), P, Q, constant_mass=true)

  # solve with SSP RK 3
  t0, tF = 0.0, 2*π
  _dt = dx(nc(panel_model))*CFL/p_fe
  dt = floor(_dt,sigdigits=1)


  solver = RungeKutta(LUSolver(), LUSolver(), dt, :EXRK_SSP_3_3)
  solT = solve(solver, opT, t0, tF, uh0)

  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
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


function transient_advection_dg_errors(panel_model,args...)
  e_u  = transient_advection_dg(panel_model,args...)
  return e_u,false,false
end

function transient_advection_dg_convergence_test(n_ref_lvls,u,vX,CFL=0.1,return_vtk=false)

  for p_fe in [1]
    errs,ns,dxs,slope = convergence_test(transient_advection_dg_errors,n_ref_lvls,u,vX,p_fe,CFL,return_vtk)
    plot()
    plot_convergence(errs,ns,dxs,slope;
        leginf=["u: p=$p_fe"],
        colors=[palette(:tab10)[p_fe]],
        ls=[:solid, :dot], )
  end
  savefig(plotsdir()*"/transient_advection_dg_convergence")

end
