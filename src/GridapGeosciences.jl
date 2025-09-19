module GridapGeosciences
using DrWatson
using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using Gridap.Algebra, Gridap.FESpaces
using LinearAlgebra
using Plots, LaTeXStrings
using FillArrays
using Test
using JLD2

dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)
!isdir(plotsdir()) && mkdir(plotsdir())

export dir

const global RADIUS = 1.0
export RADIUS

include("Helpers/Helpers.jl")

include("Fields/Fields.jl")

include("Geometry/Geometry.jl")

include("Adaptivity/Adaptivity.jl")

include("ODEs/ODEs.jl")

include("Visualisation/Visualisation.jl")


# include("Lauritzen_functions.jl")
# export gaussian_hill, cosine_bell,slotted_cylinders, correlated_cosine_bell
# export nondivergent_velocity, divergent_velocity

# include("Williamson_functions_v2.jl")
# export xyz2θϕr, θϕ2xyz, spherical_to_cartesian_matrix
# export f₀, u₀, h₀, η₀, q₀, topography, _topography

include("Exports.jl")

end
