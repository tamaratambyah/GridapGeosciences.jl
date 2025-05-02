using DrWatson
using Gridap


include("../src/initialise.jl")
coarse_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
latlon_model = LatlonDiscreteModel(coarse_model)




Ω = Triangulation(coarse_model)
cf = CellField(x->x[1],Ω)
pts = get_cell_points(Ω)

Ω_latlon = Triangulation(latlon_model)
cmap_latlon = map(x-> γinv ∘ (σ ∘ (PanelRotationField(rp1[x]) ∘ σ) ), 1:6) # latlon_p -> ambient -> parametric
cf_mapped = lazy_map((∘),get_data(cf),cmap_latlon)
cf_latlon = CellData.GenericCellField(cf_mapped,Ω_latlon,PhysicalDomain() )



cell_vals = cf_latlon(get_cell_points(Ω_latlon))
