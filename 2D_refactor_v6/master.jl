using DrWatson
using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers
using LinearAlgebra

dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

global RADIUS = 1.0

include("CellData/CellFields.jl")

include("Fields/MatMultField.jl")
include("Fields/ForwardMap.jl")
include("Fields/InverseMap.jl")

include("Geometry/cube_surface.jl")
include("Geometry/panel_matrices.jl")
include("Geometry/parametric_model.jl")
include("Geometry/ambient_model.jl")

include("Visualisation/writevtk.jl")

include("panel_ids_from_refinement.jl")
include("analytical_functions.jl")


cube_model = coarse_cube_model(π/4,6)
cube_model = Gridap.Adaptivity.refine(cube_model)
writevtk(Triangulation(cube_model),dir*"/cube_model",append=false)

panel_model,panel_ids = parametric_model(cube_model)

sphere_model = ambient_model(panel_model,panel_ids)
writevtk(Triangulation(sphere_model),dir*"/ambient_model",append=false)
