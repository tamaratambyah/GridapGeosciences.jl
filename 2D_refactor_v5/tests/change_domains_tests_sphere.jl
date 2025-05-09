using Gridap
include("../src/initialise.jl")

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
ambient_model = get_ambient_model(manifold_model)
panel_ids = get_panel_ids(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)

pts_parametric = get_cell_points(Ω_parametric)
pts_ambient = get_cell_points(Ω_ambient)

################################################################################
#### functions in the parametric space
################################################################################
cmap_ambient = map(x-> InvGnomonicField() ∘ SigmaField(r) ∘ PanelRotationField(rp1_3D[x]), panel_ids)  # ambient -> parametric


alpha(x) = x[1]

cf_parametric = CellField(alpha,Ω_parametric)
cvals_parametric = cf_parametric(pts_parametric)

cf_mapped = lazy_map(Broadcasting(∘),get_data(cf_parametric),cmap_ambient)
cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_ambient = cf_ambient(get_cell_points(Ω_ambient))

display(cvals_parametric)
display(cvals_ambient)
cvals_parametric ≈ cvals_ambient

#####
beta(x) = x[2]

cf_parametric = CellField(beta,Ω_parametric)
cvals_parametric = cf_parametric(pts_parametric)

cf_mapped = lazy_map(Broadcasting(∘),get_data(cf_parametric),cmap_ambient)
cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )
cvals_ambient = cf_ambient(get_cell_points(Ω_ambient))

display(cvals_parametric)
display(cvals_ambient)
cvals_parametric ≈ cvals_ambient



writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)
writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)
