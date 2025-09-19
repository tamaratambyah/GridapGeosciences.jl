module Geometry
using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using Gridap.Algebra, Gridap.FESpaces
using LinearAlgebra
using FillArrays

# using GridapGeosciences.Adaptivity
# import GridapGeosciences.Adaptivity: get_panel_ids

using GridapGeosciences.Fields
import GridapGeosciences.Fields: MatMultField, forward_jacobian

using GridapGeosciences.Helpers
import GridapGeosciences.Helpers: analytic_inv_metric, _analytic_inv_metric

include("panel_ids.jl")
include("cube_surface.jl")
include("panel_matrices.jl")
include("parametric_model.jl")
include("ambient_model.jl")
include("BoundaryTriangulations.jl")
include("SkeletonTriangulations.jl")
include("AdaptedTriangulations.jl")
include("panelwise_cellfield.jl")

export get_panel_ids
export pullback_area_form
export pushforward_normal, get_facet_normal, get_mapped_facet_normal
export BoundaryTriangulation
export generate_ptr, coarse_cube_model

export coarse_parametric_model
export R1p, A_cube2panel, A_panel2cube, b_panel2cube
export ParametricDiscreteModel

export get_nodes_from_coords

export panelwise_cellfield



end
