module Geometry
using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using Gridap.Algebra, Gridap.FESpaces
using LinearAlgebra
using FillArrays

# using GridapGeosciences.Adaptivity
# import GridapGeosciences.Adaptivity: get_panel_ids

import Gridap.Geometry: TriangulationView

using GridapGeosciences.Fields
import GridapGeosciences.Fields: MatMultField

using GridapGeosciences.Helpers
import GridapGeosciences.Helpers: inv_metric, forward_jacobian

include("panel_ids.jl")
include("cube_surface.jl")
include("panel_matrices.jl")
include("parametric_model.jl")
include("ambient_model.jl")
include("BoundaryTriangulations.jl")
include("SkeletonTriangulations.jl")
include("AdaptedTriangulations.jl")
include("ParametricCellField.jl")
include("TriangulationView.jl")
include("TriangulationPanelIds.jl")

export get_panel_ids, get_forward_map_generator, geo_map_func, latlon_geo_map_func
export pullback_area_form
export pushforward_normal, get_facet_normal, get_mapped_facet_normal
export BoundaryTriangulation
export generate_ptr, coarse_cube_model

export coarse_parametric_model
export R1p, A_cube2panel, A_panel2cube, b_panel2cube
export ParametricDiscreteModel

export get_nodes_from_coords

export ParametricCellField
export _pushforward_normal
export _pullback_area_form

export get_radius, get_thickness

end
