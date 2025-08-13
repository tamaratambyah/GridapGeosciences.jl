using DrWatson
using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers
using LinearAlgebra
using Plots, LaTeXStrings

dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)
!isdir(plotsdir()) && mkdir(plotsdir())


global RADIUS = 1.0

include("Adaptivity/panel_ids_from_refinement.jl")

include("CellData/CellFields.jl")

include("Fields/ForwardMap.jl")
include("Fields/InverseMap.jl")
include("Fields/MatMultField.jl")

include("Geometry/cube_surface.jl")
include("Geometry/panel_matrices.jl")
include("Geometry/parametric_model.jl")
include("Geometry/ambient_model.jl")

include("Helpers/helpers.jl")
include("Helpers/convergence_tools.jl")
include("Helpers/overloads.jl")

include("Visualisation/writevtk.jl")

include("analytical_functions.jl")
include("coordinate_mappings.jl")


cube_model = coarse_cube_model(π/4,6)
cube_model = Gridap.Adaptivity.refine(cube_model)
writevtk(Triangulation(cube_model),dir*"/cube_model",append=false)

panel_model,panel_ids = parametric_model(cube_model)

sphere_model = ambient_model(panel_model,panel_ids)
writevtk(Triangulation(sphere_model),dir*"/ambient_model",append=false)



panel_model = parametric_model(cube_model)
sphere_model = ambient_model(panel_model)
