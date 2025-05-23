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

function f(X)
  VectorValue(-X[2],X[1],0)
end


writevtk(Ω_ambient,dir*"/ambient_tangent",cellfields=["u"=>f],append=false)
