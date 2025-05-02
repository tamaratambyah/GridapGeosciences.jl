using DrWatson
using Gridap
using GridapMakie, GLMakie, FileIO


include("../src/initialise.jl")
coarse_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)

parametric_grid = get_grid(coarse_model)
ambient_grid = get_ambient_grid(parametric_grid)
_grid = UnstructuredGrid(
  get_node_coordinates(ambient_grid),
  get_cell_node_ids(parametric_grid),
  get_reffes(parametric_grid),
  get_cell_type(parametric_grid),OrientationStyle(parametric_grid),
  nothing,
  get_cell_map(parametric_grid)
)

_model = UnstructuredDiscreteModel(_grid,get_grid_topology(coarse_model),get_face_labeling(coarse_model)) |>simplexify

latlon_model = LatlonDiscreteModel(coarse_model) |>simplexify

panel_ids = get_panel_ids(latlon_model)
Ω = Triangulation(_model)
cf = CellField(x->x[1],Ω)

_Ω = Triangulation(latlon_model)

cmap_latlon = map(x-> γinv ∘ (σ ∘ (PanelRotationField(rp1[x]) ∘ σ) ), [1:6...,1:6...]) # latlon_p -> ambient -> parametric
cf_mapped = lazy_map((∘),get_data(cf),cmap_latlon)
cf_latlon = CellData.GenericCellField(cf_mapped,_Ω,PhysicalDomain() )



fig, _ , plt = plot(_Ω, cf_latlon)
wireframe!(_Ω, color=:black, linewidth=2)
