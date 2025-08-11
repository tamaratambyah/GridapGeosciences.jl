function cartesian_to_latlon(XYZ)
  X,Y,Z = XYZ
  θ = atan(Y,X)
  ϕ = atan(Z,sqrt(X^2 + Y^2))
  Point(θ,ϕ)
end

function latlon_to_cartesian(θϕ)
  θ,ϕ = θϕ
  X = cos(θ)*cos(ϕ)
  Y = sin(θ)*cos(ϕ)
  Z = sin(ϕ)
  RADIUS*Point(X,Y,Z)
end

function latlon_to_panel(p)
  function _latlon_to_panel(θϕ)
    XYZ = latlon_to_cartesian(θϕ)
    αβ = inverse_map(XYZ,p)
    return αβ
  end
end


function panel_to_latlon(αβ,p)
  XYZ = forward_map(αβ,p)
  θϕ = cartesian_to_latlon(XYZ)
  return θϕ
end


uθϕ(θϕ) = sin(θϕ[2])
uθϕ2(θϕ) = cos(θϕ[2])
uX(XYZ) =  XYZ[3]/RADIUS
uX2(XYZ) = cos( XYZ[3]/RADIUS )

function uX_panel(p,u)
  function _u(αβ)
    XYZ = forward_map(αβ,p)
    u(XYZ)
  end
end



function uθϕ_panel(p,u)
  function _u(αβ)
    θϕ = panel_to_latlon(αβ, p)
    u(θϕ)
  end
end

Ω_panel = Triangulation(panel_model)
Ω_sphere = Triangulation(ambient_model)

# ################################################################################
# ### scalar function
# ################################################################################
# ### on panel
# cell_field = map(p->GenericField(uX_panel(p,uX)),panel_ids)
# uX_cf_panel =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

# cell_field = map(p->GenericField(uθϕ_panel(p,uθϕ)),panel_ids)
# uθϕ_cf_panel =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())


# #### on sphere
# inv_f = lazy_map(p->InverseMap(p),panel_ids)
# ucf_mapped = lazy_map(Broadcasting(∘),get_data(uX_cf_panel),inv_f)
# ucf_ambient = CellData.GenericCellField(ucf_mapped,Triangulation(ambient_model),PhysicalDomain() )


# cell_field = map(p->GenericField(uX),panel_ids)
# ucf_sphere =  CellData.GenericCellField(cell_field,Ω_sphere,PhysicalDomain())
# writevtk(Triangulation(ambient_model),dir*"/ambient_model",
#           cellfields=["u"=>ucf_sphere,"u_mapped"=>ucf_ambient],append=false)


# ################################################################################
# ### take gradient
# ################################################################################
# ### on panel

# gradu_covarient_coeffs = gradient(uX_cf_panel)
# _gradu_covarient_coeffs = change_domain(gradu_covarient_coeffs,PhysicalDomain(),ReferenceDomain())

# cell_field = map(p->GenericField(inverse_jacobian(p)),panel_ids)
# basis_vectors =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

# gradu_phys =  basis_vectors ⋅ gradu_covarient_coeffs
# _gradu_phys = change_domain(gradu_phys,PhysicalDomain(),ReferenceDomain())

# cell_field = map(p->GenericField(uX_panel(p,uX2)),panel_ids)
# ucos =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

# function magnitude(x)
#   x[1]*x[1] + x[2]*x[2] + x[3]*x[3]
# end

# for p in collect(1:3)
#   mask = panel_ids.==p
#   Ωp = Triangulation(panel_model,mask)
#   writevtk(Ωp,dir*"/panel$(p)_model_cos",
#           cellfields=["u"=>uX_cf_panel,"gradX"=>gradient(uX_cf_panel),"gradtheta"=>gradient(uθϕ_cf_panel),
#           "basis"=>basis_vectors,"sgrad"=>gradu_phys,"cos"=>ucos,"e"=>(Operation(magnitude)(gradu_phys))-ucos],append=false)
# end



################################################################################
### interpolate into FE space, then take covariant gradient
################################################################################
### on panel
Ω_panel = Triangulation(panel_model)
Ω_sphere = Triangulation(ambient_model)

cell_field = map(p->GenericField(inverse_jacobian(p)),panel_ids)
basis_vectors =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

cell_field = map(p->GenericField(measure(p)),panel_ids)
meas_cf =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

cell_field = map(p->GenericField(inv_metric(p)),panel_ids)
inv_metric_cf = CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

cell_field = map(p->GenericField(uθϕ_panel(p,uθϕ)),panel_ids)
usin_cf_panel =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

cell_field = map(p->GenericField(uθϕ_panel(p,uθϕ2)),panel_ids)
ucos_cf_panel =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,2); conformity=:H1)
U = TrialFESpace(V)

uh_sin = interpolate(usin_cf_panel,U)
uh_cos = interpolate(ucos_cf_panel,U)

grad_sinu_phys =  basis_vectors ⋅ gradient(uh_sin)
grad_sinu_phys2 =  basis_vectors ⋅ inv_metric_cf ⋅ gradient(uh_sin)
grad_cosu_phys = basis_vectors ⋅ gradient(uh_cos)

for p in collect(1:6)
  mask = panel_ids.==p
  Ωp = Triangulation(panel_model,mask)
  writevtk(Ωp,dir*"/panel$(p)_model",
          cellfields=["uh_sin"=>uh_sin,"uh_cos"=>uh_cos,
          "grad_sinu"=>gradient(uh_sin),"grad_cosu"=>gradient(uh_cos),
          "sgrad_sin"=>grad_sinu_phys,"sgrad_cos"=>grad_cosu_phys,
          "basis"=>basis_vectors,"inv_metric"=>inv_metric_cf,"meas"=>meas_cf,
          "symmetric"=>grad_sinu_phys2
          ],append=false)
end

### map to sphere
inv_f = lazy_map(p->InverseMap(p),panel_ids)

_uh_sin = change_domain(uh_sin,ReferenceDomain(),PhysicalDomain())
cf_mapped = lazy_map(Broadcasting(∘),get_data(_uh_sin),inv_f)
uh_sin_ambient = CellData.GenericCellField(cf_mapped,Ω_sphere,PhysicalDomain() ) # ambient cell field

_uh_cos = change_domain(uh_cos,ReferenceDomain(),PhysicalDomain())
cf_mapped = lazy_map(Broadcasting(∘),get_data(_uh_cos),inv_f)
uh_cos_ambient = CellData.GenericCellField(cf_mapped,Ω_sphere,PhysicalDomain() ) # ambient cell field


_grad_sinu_phys = change_domain(grad_sinu_phys,ReferenceDomain(),PhysicalDomain())
cf_mapped = lazy_map(Broadcasting(∘),get_data(_grad_sinu_phys),inv_f)
grad_sinu_ambient = CellData.GenericCellField(cf_mapped,Ω_sphere,PhysicalDomain() ) # ambient cell field

_grad_cosu_phys = change_domain(grad_cosu_phys,ReferenceDomain(),PhysicalDomain())
cf_mapped_cos = lazy_map(Broadcasting(∘),get_data(_grad_cosu_phys),inv_f)
grad_cosu_ambient = CellData.GenericCellField(cf_mapped_cos,Ω_sphere,PhysicalDomain() )



writevtk(Triangulation(ambient_model),dir*"/ambient_model",
          cellfields=["sinu"=>uh_sin_ambient,"cosu"=>uh_cos_ambient,
          "sinuX"=>uX,"cosuX"=>uX2,
          "gradsin"=>grad_sinu_ambient,"gradcos"=>grad_cosu_ambient
          ],append=false)
