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

  # α,β = αβ

  # if p == 1
  #   θ = α
  #   ϕ = atan( tan(β)*cos(α)  )
  # elseif p == 2
  #   θ = atan( tan(α) , tan(β) )
  #   ϕ = atan( sin(atan( tan(α) , tan(β) ))  ,tan(α))
  # elseif p == 3
  #   θ = -atan(1,tan(α))
  #   ϕ = atan( tan(β)* sin(-atan(1,tan(α))) )
  # end

  # Point(θ,ϕ)

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

################################################################################
### scalar function
################################################################################
### on panel
cell_field = map(p->GenericField(uX_panel(p,uX)),panel_ids)
uX_cf_panel =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())



cell_field = map(p->GenericField(uθϕ_panel(p,uθϕ)),panel_ids)
uθϕ_cf_panel =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())


#### on sphere
inv_f = lazy_map(p->InverseMap(p),panel_ids)
ucf_mapped = lazy_map(Broadcasting(∘),get_data(uX_cf_panel),inv_f)
ucf_ambient = CellData.GenericCellField(ucf_mapped,Triangulation(ambient_model),PhysicalDomain() )


cell_field = map(p->GenericField(uX),panel_ids)
ucf_sphere =  CellData.GenericCellField(cell_field,Ω_sphere,PhysicalDomain())
writevtk(Triangulation(ambient_model),dir*"/ambient_model",
          cellfields=["u"=>ucf_sphere,"u_mapped"=>ucf_ambient],append=false)


################################################################################
### take gradient
################################################################################
### on panel

gradu_covarient_coeffs = gradient(uX_cf_panel)
_gradu_covarient_coeffs = change_domain(gradu_covarient_coeffs,PhysicalDomain(),ReferenceDomain())

cell_field = map(p->GenericField(inverse_jacobian(p)),panel_ids)
basis_vectors =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

gradu_phys =  basis_vectors ⋅ gradu_covarient_coeffs
_gradu_phys = change_domain(gradu_phys,PhysicalDomain(),ReferenceDomain())

cell_field = map(p->GenericField(uX_panel(p,uX2)),panel_ids)
ucos =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

function magnitude(x)
  x[1]*x[1] + x[2]*x[2] + x[3]*x[3]
end

for p in collect(1:3)
  mask = panel_ids.==p
  Ωp = Triangulation(panel_model,mask)
  writevtk(Ωp,dir*"/panel$(p)_model_cos",
          cellfields=["u"=>uX_cf_panel,"gradX"=>gradient(uX_cf_panel),"gradtheta"=>gradient(uθϕ_cf_panel),
          "basis"=>basis_vectors,"sgrad"=>gradu_phys,"cos"=>ucos,"e"=>(Operation(magnitude)(gradu_phys))-ucos],append=false)
end



################################################################################
### interpolate into FE space, then
### take gradient
################################################################################
### on panel

# cell_field = map(p->GenericField(uX_panel(p,uX)),panel_ids)
cell_field = map(p->GenericField(uθϕ_panel(p,uθϕ)),panel_ids)
usin_cf_panel =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

cell_field = map(p->GenericField(uθϕ_panel(p,uθϕ2)),panel_ids)
ucos_cf_panel =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,2); conformity=:H1)
U = TrialFESpace(V)

uh = interpolate(usin_cf_panel,U)

gradu_phys =  basis_vectors ⋅ gradient(uh)


for p in collect(1:3)
  mask = panel_ids.==p
  Ωp = Triangulation(panel_model,mask)
  writevtk(Ωp,dir*"/panel$(p)_model",
          cellfields=["uh"=>uh,"gradu"=>gradient(uh),"sgrad"=>gradu_phys,"cos"=>ucos_cf_panel],append=false)
end

_gradu_phys = change_domain(gradu_phys,ReferenceDomain(),PhysicalDomain())

cf_mapped = lazy_map(Broadcasting(∘),get_data(_gradu_phys),inv_f)
gradu_phys_ambient = CellData.GenericCellField(cf_mapped,Ω_sphere,PhysicalDomain() ) # ambient cell field

gradcos = basis_vectors ⋅ gradient(ucos)
cf_mapped_cos = lazy_map(Broadcasting(∘),get_data(gradcos),inv_f)
gradcos_ambient = CellData.GenericCellField(cf_mapped_cos,Ω_sphere,PhysicalDomain() )

writevtk(Triangulation(ambient_model),dir*"/ambient_model",
          cellfields=["gradu"=>gradu_phys_ambient],append=false)

writevtk(Triangulation(ambient_model),dir*"/ambient_model",
          cellfields=["u"=>ucf_sphere,"u_mapped"=>ucf_ambient,"gradu"=>gradu_phys_ambient,"gradcos"=>gradcos_ambient],append=false)
