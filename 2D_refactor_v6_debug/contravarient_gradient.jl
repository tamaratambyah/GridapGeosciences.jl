################################################################################
### interpolate into FE space, then take covariant gradient
################################################################################
### on panel
Ω_panel = Triangulation(panel_model)
Ω_sphere = Triangulation(ambient_model)

cell_field = map(p->GenericField(forward_jacobian(p)),panel_ids)
basis_vectors =  CellData.GenericCellField(cell_field,Ω_panel,PhysicalDomain())

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

grad_sinu_phys =  basis_vectors ⋅ (inv_metric_cf ⋅ gradient(uh_sin))
grad_cosu_phys = basis_vectors ⋅ (inv_metric_cf ⋅ gradient(uh_cos))

for p in collect(1:6)
  mask = panel_ids.==p
  Ωp = Triangulation(panel_model,mask)
  writevtk(Ωp,dir*"/panel$(p)_model_contra",
          cellfields=["uh_sin"=>uh_sin,"uh_cos"=>uh_cos,
          "grad_sinu"=>gradient(uh_sin),"grad_cosu"=>gradient(uh_cos),
          "sgrad_sin"=>grad_sinu_phys,"sgrad_cos"=>grad_cosu_phys,
          "basis"=>basis_vectors,"inv_metric"=>inv_metric_cf
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



writevtk(Triangulation(ambient_model),dir*"/ambient_model_contr",
          cellfields=["sinu"=>uh_sin_ambient,"cosu"=>uh_cos_ambient,
          "gradsin"=>grad_sinu_ambient,"gradcos"=>grad_cosu_ambient,
          "esin"=>grad_sinu_ambient-uh_cos_ambient
          ],append=false)
