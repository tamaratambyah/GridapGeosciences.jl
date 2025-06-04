"""
u_scalar_latlon2parametric

Map analytical scalar valued function from θ,ϕ -> α,β
1. Map a point α,β -> θ,ϕ
2. Evaluate g(α,β) = f(θ,ϕ)
"""
function u_scalar_latlon2parametric(p::Int,uθϕ::Function)
  function _u(αβ)
    latlon_panel1 = evaluate(GnomonicMap(), αβ)
    sphere_panel1 = evaluate(Sigma(),latlon_panel1)
    sphere_panelp = evaluate(PanelRotationMap(r1p_3D[p]), sphere_panel1)
    latlon_panelp = evaluate(Sigma(),sphere_panelp)
    uθϕ(latlon_panelp)
  end
end


biform(u,v,dΩ) = ∫(v*u)dΩ
liform(v,f,dΩ)   = ∫(v*f)dΩ

function parametric_cf_2_ambient(manifold_model,degree::Int,cf_parametric::CellField)
  ambient_model = get_ambient_model(manifold_model)
  panel_ids = get_panel_ids(manifold_model)

  Ω_parametric = Triangulation(manifold_model)
  Ω_ambient = Triangulation(ambient_model)

  inv_mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)

  # evaluate at quadrature points
  uh_parametric = L2_scalar_cf_projection(manifold_model,degree,cf_parametric)

  ## map parametric FEFunction back to ambient space
  _uh = change_domain(uh_parametric,ReferenceDomain(),PhysicalDomain())

  cf_mapped = lazy_map(Broadcasting(∘),get_data(_uh),inv_mapping)
  cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )

  uh_ambient = L2_scalar_cf_projection(ambient_model,degree,cf_ambient)

  return uh_parametric, uh_ambient
end


function L2_scalar_cf_projection(model::DiscreteModel,degree::Int,cf::CellField)
  Ω = Triangulation(model)

  # evaluation at quadrature points
  dΩ = Measure(Ω,degree)
  L2 = FESpace(model,ReferenceFE(lagrangian,Float64,1), conformity=:L2)
  op = AffineFEOperator((u,v)->biform(u,v,dΩ),
                        v->liform(v,cf,dΩ),
                        L2,L2)
  uh = solve(LUSolver(),op)
  return uh
end
