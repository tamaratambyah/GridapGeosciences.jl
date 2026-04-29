"""
test mass conservation by solving
u + ∇ᵧ(φ) = u₀
∇ᵧ⋅u = 0
for arbitray vector field u₀
"""


function mass_conservation(panel_model,p_fe::Int,dir::String,func::Function,scalar_field::Bool,return_vtk=false)
  lvl = nref(nc(panel_model))
  println("p_fe = $(p_fe); nref = $lvl")

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  vec_contra_cf = if scalar_field
    ParametricCellField(contr_gradf(func),Ω_panel,panel_ids)
  else
    ParametricCellField(contra_v(func),Ω_panel,panel_ids)
  end

  metric_cf = ParametricCellField(metric,Ω_panel,panel_ids)
  meas_cf = ParametricCellField(sqrtg,Ω_panel,panel_ids)
  grad_meas_cf = ParametricCellField(grad_meas,Ω_panel,panel_ids)
  covariant_basis_cf = ParametricCellField(covariant_basis,Ω_panel,panel_ids)

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
    cell_geo_map = geo_map_func(Ω_panel)

    u_proj = covariant_basis_cf ⋅ vec_contra_cf
    u_projh = covariant_basis_cf ⋅ uh

    panel_cfs = [ u_proj, u_projh, ph]
    labels = ["u0", "u_projh", "p"]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  println("initial divergence: $s_div0")

  return abs(s_div), false,false

end



function mass_conservation_convergence_test(ranks::AbstractArray,nprocs::Int,analytic_funcs,scalar_field::Bool,
  n_ref_lvls=4,ps=[1],return_vtk=false)
  # serial models
  models  = get_refined_models(n_ref_lvls)
  dir = datadir("MassConservationConvergence")
  !isdir(dir) && mkdir(dir)

  for (key, val) in analytic_funcs
    _dir = dir*"/func_$(key)"
    !isdir(_dir) && mkdir(_dir)

    errors = Vector{Vector{Float64}}(undef,length(ps))
    ns = Vector{Vector{Float64}}(undef,length(ps))
    dxs = Vector{Vector{Float64}}(undef,length(ps))
    slopes = Vector{Float64}(undef,length(ps))
    for (i,p_fe) in enumerate(ps)
      errors[i],ns[i],dxs[i],slopes[i] = h_convergence_test(models,mass_conservation,p_fe,_dir,val,scalar_field,return_vtk)
    end

    @test all(all.(map(x->x.< 1e-12,errors)))

    output = @strdict errors ns dxs slopes ps
    safesave(datadir(_dir, ("convergence.jld2")), output)
    plot_convergence_from_saved(_dir,"convergence",["p"])
  end

end
