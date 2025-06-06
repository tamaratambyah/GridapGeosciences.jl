using Gridap
include("../src/initialise.jl")


manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(manifold_model)

ambient_model = get_ambient_model(manifold_model)
latlon_model = get_latlon_model(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)
Ω_latlon = Triangulation(latlon_model)

pts_parametric = get_cell_points(Ω_parametric)
pts_ambient = get_cell_points(Ω_ambient)
pts_latlon = get_cell_points(Ω_latlon)


################################################################################
##### analytical function in ambient space
################################################################################

RT_ambient = FESpace(ambient_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)
analytic_u_ambient = interpolate(tangent_f,RT_ambient)

# writevtk(Ω_ambient,dir*"/ambient_tangent",cellfields=["u"=>tangent_f,"uh"=>ambient_uh],append=false)


################################################################################
##### Interoplate vector valued function
################################################################################

RT = FESpace(manifold_model,ReferenceFE(raviart_thomas,Float64,1), conformity=:HDiv)

## interpolate everywhere
free_values = zero_free_values(RT)
dirichlet_values = zero_dirichlet_values(RT)

## get cell vals
s = get_fe_dof_basis(RT)
trian = get_triangulation(s)
cell_field = map(p->GenericField(u_vector_ambient2parametric(p,tangent_f)),panel_ids)
f = CellData.GenericCellField(cell_field,trian,PhysicalDomain())
cell_vals = s(f)

## interpolate!
gather_free_values!(free_values,RT,cell_vals)
uh = FEFunction(RT,free_values)
# writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>uh],append=false)

################################################################################
##### Map back to ambient space for visualisation
################################################################################
mapping = map(x-> PanelRotationField(r1p_3D[x]) ∘ SigmaField(RADIUS) ∘ GnomonicField() , panel_ids)
inv_mapping = map(x-> InvGnomonicField() ∘ InvSigmaField(RADIUS) ∘ PanelRotationField(rp1_3D[x]), panel_ids)


Jt = lazy_map(Broadcasting(gradient),mapping)
J = lazy_map(Operation(transpose),Jt)

_uh = change_domain(uh.cell_field,ReferenceDomain(),PhysicalDomain())

_cf_mapped = lazy_map(Broadcasting(⋅),J,get_data(_uh))
cf_mapped = lazy_map(Broadcasting(∘),_cf_mapped,inv_mapping)

cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_ambient = cf_ambient(pts_ambient)

### Convert the cellfield into a FEFunction by evaluating the field at the dofs
############ get cell vals on ambient RT
ambient_free_values = zero_free_values(RT_ambient)
ambient_s = get_fe_dof_basis(RT_ambient)
ambient_cell_vals = ambient_s(cf_ambient)

## interpolate!
gather_free_values!(ambient_free_values,RT_ambient,ambient_cell_vals)
uh_ambient = FEFunction(RT_ambient,ambient_free_values)

writevtk(Ω_ambient,dir*"/ambient_tangent",
      cellfields=["f"=>tangent_f,"u_ambient_cf"=>cf_ambient,"uh_ambient"=>uh_ambient,
      "eh"=>analytic_u_ambient-uh_ambient,
      "u_analytic_ambient"=>analytic_u_ambient],append=false)


e = analytic_u_ambient - uh_ambient
dΩ = Measure(Ω_ambient,2)
l2(e,dΩ)
