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

const global RADIUS = 1.0
export RADIUS

include("Helpers/Helpers.jl")

include("Fields/Fields.jl")

include("Geometry/Geometry.jl")

include("Adaptivity/Adaptivity.jl")

include("ODEs/ODEs.jl")

include("Visualisation/Visualisation.jl")

include("convergence_tools.jl")
export l2
export get_refined_models, convergence_test, h_convergence_test
export plot_convergence, plot_error, convergence_rate, plot_convergence_from_saved
export print_convergence_results
export nc, dx, nref

include("Distributed/Distributed.jl")

include("Exports.jl")

end
