using Gridap
include("../src/initialise.jl")

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(1),cube)
ambient_model = get_ambient_model(manifold_model)
panel_ids = get_panel_ids(manifold_model)

Ω_parametric = Triangulation(manifold_model)
Ω_ambient = Triangulation(ambient_model)

pts_parametric = get_cell_points(Ω_parametric)
pts_ambient = get_cell_points(Ω_ambient)

################################################################################
#### functions in the parametric space
################################################################################
alpha(x) = x[1]
beta(x) = x[2]

cmap_ambient = map(x-> BumpField( A_bump,B_bump,b_bump) ∘ PanelRotationField(rp1_3D[x]), panel_ids)  # ambient -> parametric


cf_parametric = CellField(beta,Ω_parametric)
cf_mapped = lazy_map(Broadcasting(∘),get_data(cf_parametric),cmap_ambient)
cf_ambient = CellData.GenericCellField(cf_mapped,Ω_ambient,PhysicalDomain() )

cvals_parametric = cf_parametric(pts_parametric)
cvals_ambient = cf_ambient(pts_ambient)

cvals_parametric == cvals_ambient

writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)
writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)


################################################################################
#### functions in the ambient space
################################################################################
f(x) = (x[1]+x[2]+x[3])^2

cmap_parametric = map(x-> PanelRotationField(r1p_3D[x]) ∘ BumpField( A_bump,B_bump,b_bump), panel_ids)  # parametric -> ambient

cf_ambient = CellField(f,Ω_ambient)
cf_mapped = lazy_map(Broadcasting(∘),get_data(cf_ambient),cmap_parametric)
cf_parametric = CellData.GenericCellField(cf_mapped,Ω_parametric,PhysicalDomain() )

cvals_parametric = cf_parametric(pts_parametric)
cvals_ambient = cf_ambient(pts_ambient)

cvals_parametric == cvals_ambient


writevtk(Ω_parametric,dir*"/parametric",cellfields=["u"=>cf_parametric],append=false)
writevtk(Ω_ambient,dir*"/ambient",cellfields=["u"=>cf_ambient],append=false)
