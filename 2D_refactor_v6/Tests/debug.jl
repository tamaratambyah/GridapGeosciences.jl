
####################### vector mapping
cube_model = coarse_cube_model(π/4,6)
cube_model = Gridap.Adaptivity.refine(cube_model)
cube_model = Gridap.Adaptivity.refine(cube_model)

panel_model = parametric_model(cube_model)
# panel_model = coarse_parametric_model()
# panel_model = Gridap.Adaptivity.refine(panel_model)
# panel_model = Gridap.Adaptivity.refine(panel_model)

Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,20)
panel_ids = get_panel_ids(panel_model)

vX(XYZ) = VectorValue(-1.0*XYZ[2],XYZ[1],0.0)
vX(XYZ) = VectorValue(XYZ[1]*XYZ[2],XYZ[2]*XYZ[3],XYZ[3]^2-RADIUS^2)

_vX(p) = αβ -> vX(forward_map(p)(αβ))

# forward_pinv_jacobian(p) = αβ -> forward_pinv_jacobian(αβ,p)




function contra_v(αβ,p)

  α,β = αβ
  rho = sqrt(1 + (tan(α))^2 + (tan(β))^2 )
  drho_da = - tan(α)*(sec(α))^2 / ( rho^3 )
  drho_db = - tan(β)*(sec(β))^2 / ( rho^3 )

  if p == 1
    X = 1/rho
    Y = 1/rho * tan(α)
    Z = 1/rho * tan(β)

    dXda = drho_da
    dXdb = drho_db
    dYda = drho_da*tan(α) + 1/rho*(sec(α))^2
    dYdb = drho_db*tan(α)
    dZda = drho_da*tan(β)
    dZdb = drho_db*tan(β) + 1/rho*(sec(β))^2
  elseif p == 2
    X = -1/rho * tan(β)
    Y = 1/rho * tan(α)
    Z = 1/rho

    dXda = -( drho_da*tan(β) )
    dXdb = -( drho_db*tan(β) + 1/rho*(sec(β))^2 )
    dYda = drho_da*tan(α) + 1/rho*(sec(α))^2
    dYdb = drho_db*tan(α)
    dZda = drho_da
    dZdb = drho_db
  elseif p == 3
    X = -1/rho * tan(α)
    Y = 1/rho
    Z = 1/rho * tan(β)

    dXda = -( drho_da*tan(α) + 1/rho*(sec(α))^2  )
    dXdb = -( drho_db*tan(α) )
    dYda = drho_da
    dYdb = drho_db
    dZda = drho_da*tan(β)
    dZdb = drho_db*tan(β) + 1/rho*(sec(β))^2
  elseif p == 4
    X = -1/rho
    Y = 1/rho * tan(β)
    Z = 1/rho * tan(α)

    dXda = -drho_da
    dXdb = -drho_db
    dYda = drho_da*tan(β)
    dYdb = drho_db*tan(β) + 1/rho*(sec(β))^2
    dZda = drho_da*tan(α) + 1/rho*(sec(α))^2
    dZdb = drho_db*tan(α)
  elseif p == 5
    X = -1/rho * tan(α)
    Y = 1/rho * tan(β)
    Z = -1/rho

    dXda = -(drho_da*tan(α) + 1/rho*(sec(α))^2)
    dXdb = -(drho_db*tan(α) )
    dYda = drho_da*tan(β)
    dYdb = drho_db*tan(β) + 1/rho*(sec(β))^2
    dZda = -drho_da
    dZdb = -drho_db
  elseif p == 6
    X = -1/rho * tan(β)
    Y = -1/rho
    Z = 1/rho * tan(α)

    dXda = -(drho_da*tan(β))
    dXdb = -(drho_db*tan(β) + 1/rho*(sec(β))^2)
    dYda = -drho_da
    dYdb = -drho_db
    dZda = drho_da*tan(α) + 1/rho*(sec(α))^2
    dZdb = drho_db*tan(α)

  end

  ## J = [dXda dXdb
  ##      dYda dYdb
  ##      dZda dZdb  ]
  ## As a TensorValue data = (dXda,dYda,dZda, dXdb, dYdb, dZdb)

  J = RADIUS*TensorValue{3,2}(dXda,dYda,dZda, dXdb,dYdb,dZdb)
  Jt = transpose(J)#RADIUS*TensorValue{2,3}(dXda,dXdb, dYda,dYdb, dZda,dZdb)
  XYZ = RADIUS*VectorValue(X,Y,Z)
  (inv(Jt⋅J)⋅Jt) ⋅ vX(XYZ)
  # (analytic_inv_metric(αβ) ⋅ Jt )  ⋅ vX(XYZ)
end
contra_v(p) = αβ -> contra_v(αβ,p)

# contra_v(p) = αβ -> forward_pinv_jacobian(p)(αβ)⋅ _vX(p)(αβ)
# contra_v(p) = αβ -> forward_pinv_jacobian(p)(αβ)⋅ _vX(p)(αβ)
projection_vX(p) = αβ -> forward_jacobian(p)(αβ) ⋅ contra_v(p)(αβ)

vec_phys = panelwise_cellfield(_vX,Ω_panel,panel_ids)
vec_project = panelwise_cellfield(projection_vX,Ω_panel,panel_ids)

vec_contra_cf = panelwise_cellfield(contra_v,Ω_panel,panel_ids)
jacobian_cf = panelwise_cellfield(forward_jacobian,Ω_panel,panel_ids)

p_fe = 1
V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
U = TrialFESpace(V)
# biform(u,v) = ∫( u⋅v )dΩ
# liform(v) = ∫( vec_contra_cf⋅v )dΩ
# op = AffineFEOperator(biform,liform,U,V)
# vec_contra_h_solve = solve(LUSolver(),op)

vec_contra_h = interpolate(vec_contra_cf,U)
project_h = jacobian_cf ⋅vec_contra_h

# e = vec_phys-vec_project
e = project_h - vec_project
# e = vec_contra_h - vec_contra_h_solve
l2(e,dΩ)

cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)

panel_cfs = [vec_phys, vec_project,e, project_h]
labels = ["u","u_proj", "e" ,"u_projh"]

cellfields = map((x,y) -> x=>y, labels,panel_cfs)

writevtk(Ω_panel,dir*"/ambient_model",cellfields=cellfields,append=false,geo_map=cell_geo_map)










# function pinvJ(J::MultiValue{Tuple{D1,D2}}) where {D1,D2}
#   Jt = transpose(J)
#   inv(Jt⋅J)⋅Jt
# end

# vX(XYZ) = VectorValue(-1.0*XYZ[2],XYZ[1],0.0)
vX(XYZ) = VectorValue(XYZ[1]*XYZ[2],XYZ[2]*XYZ[3],XYZ[3]^2-RADIUS^2)

_vX(p) = αβ -> vX(forward_map(p)(αβ))
vec_phys = panelwise_cellfield(_vX,Ω_panel,panel_ids)


# forward_map(p) = αβ -> forward_map(αβ,p)
# v_tangent(p) = αβ ->  pinvJ(forward_jacobian(p)(αβ)) ⋅ vX(forward_map(p)(αβ))

forward_pinv_jacobian(p) = αβ -> forward_pinv_jacobian(αβ,p)
# v_tangent(p) = αβ ->  forward_pinv_jacobian(p)(αβ) ⋅ _vX(p)(αβ)
# _v_tangent(p) = αβ ->  inverse_jacobian(p)(αβ) ⋅ _vX(p)(αβ)

v_tangent(p) = αβ -> (R1p[p]⋅analytic_J1(αβ)) ⋅ ( (analytic_inv_metric(αβ) ⋅  transpose(R1p[p]⋅analytic_J1(αβ)))⋅ _vX(p)(αβ) )


# contravarient_comps = panelwise_cellfield(v_tangent,Ω_panel,panel_ids)

# covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

# vec_panel = covarient_basis_cf ⋅  contravarient_comps
vec_panel = panelwise_cellfield(v_tangent,Ω_panel,panel_ids)

e = vec_phys-vec_panel

dΩ = Measure(Ω_panel,8)
l2(e,dΩ)


cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)

panel_cfs = [vec_phys, vec_panel, e,contravarient_comps]
labels = ["u","u_panel", "e","u_con"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)

writevtk(Ω_panel,dir*"/ambient_model",cellfields=cellfields,append=false,geo_map=cell_geo_map)


writevtk_ambient(panel_model,panel_cfs,labels)
