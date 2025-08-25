
function vector_field(panel_model,vX::Function)
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)

  vec_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  return vec_contra_cf
end

function sgrad_vector_field(panel_model,f::Function)
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)

  vec_contra_cf = panelwise_cellfield(contr_gradf(f),Ω_panel,panel_ids)
  return vec_contra_cf
end

function mass_conservation(panel_model,vec_contra_cf::CellField,p_fe::Int,return_vtk=false )
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

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
    writevtk(Ω_panel,dir*"/ambient_model_nref$lvl",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  return s_div, s_div0, panel_div

end


function mass_conservation_errors(panel_model,func::Function,p_fe::Int,return_vtk=false,scalar_field=false)
  vec_contra_cf = vector_field(panel_model,func)

  if scalar_field
    vec_contra_cf = sgrad_vector_field(panel_model,func)
  end

  s_div, s_div0, = mass_conservation(panel_model,vec_contra_cf,p_fe,return_vtk )
  println("initial divergence: $s_div0")
  return abs(s_div),false
end

function mass_conservation_convergence_test(analytic_funcs,n_ref_lvls,return_vtk=false,scalar_field=false)

  for (key, val) in analytic_funcs
    plot()
    for p_fe in [1,2]
      errs,ns,dxs,slope = convergence_test(mass_conservation_errors,n_ref_lvls,val,p_fe,return_vtk,scalar_field)
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



hWilliamson(ζ,u0,ω) = θϕ -> 1 - (ω*u0 + 0.5*u0^2)*( -cos(θϕ[1])*cos(θϕ[2])*sin(ζ) +  sin(θϕ[2])*cos(ζ) )^2
ω = 1e-5
u0 = 0.1
ζ = 0.0

scalar_funcs = Dict{Symbol,Any}()
scalar_funcs[:sin] = f_sin
scalar_funcs[:XYZ] = f_XYZ
scalar_funcs[:W2] = panel_to_latlon(hWilliamson(ζ,u0,ω))

n_ref_lvls = 4
mass_conservation_convergence_test(scalar_funcs,n_ref_lvls,false,true)

#### arbitary vector fields
vecX1(XYZ) = VectorValue(XYZ[2],XYZ[3]*XYZ[2]^2,XYZ[1])

vWilliamson(ζ,u0,ω) = θϕ -> - u0*VectorValue( cos(θϕ[2])*cos(ζ) + cos(θϕ[1])*sin(θϕ[2])*sin(ζ),
                                      -sin(θϕ[1])*sin(ζ) )
vecX2 = vec_cartesian_to_latlon(vWilliamson(ζ,u0,ω))


vec_funcs = Dict{Symbol,Any}()
vec_funcs[:v] = panel_to_cartesian(tangent_vec(vecX1))
vec_funcs[:W2] = panel_to_cartesian(tangent_vec(vecX2))

n_ref_lvls = 4
mass_conservation_convergence_test(vec_funcs,n_ref_lvls,false,false)
