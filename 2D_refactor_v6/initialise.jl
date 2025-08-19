using DrWatson
using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using LinearAlgebra
using Plots, LaTeXStrings

dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)
!isdir(plotsdir()) && mkdir(plotsdir())


global RADIUS = 1.0

include("CellData/CellFields.jl")

include("Fields/ForwardMap.jl")
include("Fields/InverseMap.jl")
include("Fields/MatMultField.jl")
include("Fields/AffineField.jl")

include("Geometry/cube_surface.jl")
include("Geometry/panel_matrices.jl")
include("Geometry/parametric_model.jl")
include("Geometry/ambient_model.jl")

include("Adaptivity/panel_ids_from_refinement.jl")
include("Adaptivity/Refinement.jl")

include("Helpers/helpers.jl")
include("Helpers/convergence_tools.jl")
# include("Helpers/overloads.jl")

include("Visualisation/Vtk.jl")
include("Visualisation/helpers.jl")

include("analytical_functions.jl")
include("coordinate_mappings.jl")
include("vector_projection_analytic_functions.jl")
