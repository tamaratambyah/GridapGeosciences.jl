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


"""
u_scalar_ambient2parametric

Map analytical scalar valued function from X,Y,Z -> α,β
1. Map a point α,β -> X,Y,Z
2. Evaluate g(α,β) = f(X,Y,Z)
"""
function u_scalar_ambient2parametric(p::Int,uX::Function)
  function _u(αβ)
    # cmap = PanelRotationField(r1p_3D[p]) ∘ SigmaField(RADIUS) ∘ GnomonicField()
    # XYZ = cmap(αβ)

    θϕ = GnomonicField()(αβ)
    _XYZ = SigmaField(RADIUS)(θϕ)
    XYZ = PanelRotationField(r1p_3D[p])(_XYZ)

    uX(XYZ)
  end
end


"""
u_scalar_latlon2ambient

Map analytical scalar valued function from θ,ϕ -> X,Y,Z
1. Map a point θ,ϕ -> X,Y,Z
2. Evaluate g(θ,ϕ) = f(X,Y,Z)
"""
function u_scalar_latlon2ambient(p::Int,uθϕ::Function)
  function _u(XYZ)
    θϕ = InvSigmaField(RADIUS)(XYZ)

    uθϕ(θϕ)
  end
end


biform(u,v,dΩ) = ∫(v*u)dΩ
liform(v,f,dΩ)   = ∫(v*f)dΩ

"""
parametric_cf_2_ambient

Map a scalar cellfield defined on parametric space to a scalar cell field
defined on ambient space.

1. Change domain to physical domain
2. Apply inverse mapping so that the resulting cellfield takes points in the ambient space

Since these are scalar fields, project the resulting cellfields onto L2
"""
function parametric_cf_2_ambient(manifold_model,degree::Int,cf_parametric::CellField)
  ambient_model = get_ambient_model(manifold_model)
  panel_ids = get_panel_ids(manifold_model)

  Ω_parametric = Triangulation(manifold_model)
  Ω_ambient = Triangulation(ambient_model)

  inv_mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(RADIUS) ∘ PanelRotationField(rp1_3D[x]), panel_ids)

  # evaluate at quadrature points
  uh_parametric = L2_scalar_cf_projection(manifold_model,degree,cf_parametric)

  ## map parametric FEFunction back to ambient space
  _uh = change_domain(uh_parametric,ReferenceDomain(),PhysicalDomain())

  cf_mapped = lazy_map(Broadcasting(∘),get_data(_uh),inv_mapping)
  cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )

  uh_ambient = L2_scalar_cf_projection(ambient_model,degree,cf_ambient)

  return uh_parametric, uh_ambient
end


"""
L2_scalar_cf_projection

Computes L2 projection of a scalar cellfield, in order to evaluate at quadrature
points. This avoids any pole problem that may arise in the definition of analytic
functions.
"""
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
