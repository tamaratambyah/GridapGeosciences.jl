module  Helpers

using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using Gridap.Algebra, Gridap.FESpaces, Gridap.Helpers, Gridap.Arrays
import Gridap.TensorValues: MultiValue
using LinearAlgebra
using FillArrays
using Test # required for convergence tools (I believe it should not be here)
using GridapDistributed # required for convergence tools (I believe it should not be here)
using PartitionedArrays # required for convergence tools (I believe it should not be here)

# This using used to be in src/Helpers/forward_map_3D.jl, which is no longer in the project files
# Essentially, we need in "convergence_tools.jl" the coarse_parametric_model() function, which is defined
# in GridapGeosciences.Geometry. I still need to understand how it is possible that one can include the
# whole GridapGeoosciences package in this Helpers module, without creating circular dependencies.
using GridapGeosciences

include("ForwardMap.jl")
include("overloads.jl")
include("analytical_functions.jl")
include("analytical_functions_autodiff.jl")
include("operators.jl")
include("coordinate_mappings.jl")
include("vector_projection_analytic_functions.jl")
include("convergence_tools.jl")

export ForwardMap, ForwardMap2D, ForwardMap3D, ForwardMap2DGenerator, ForwardMap3DGenerator
export forward_jacobian, covariant_basis

export _sqrtg
export pinvJ
export sqrtg,  detg
export grad_meas
export metric, inv_metric
export perp_matrix
export surflap, surfdiv, contr_gradf, sgrad

export panel_to_cartesian, panel_to_latlon
export normal_vector_from_basis
export vector_length, normal_vec, tangent_vec
export contra_v, contra_v_comp, contra_v_perp, projection_v
export contra_v_comp3D, contra_v_perp3D
export piola
export g_star
export xyz2θϕr, θϕ2xyz, spherical_to_cartesian_matrix, xyz2θϕ

export p_convergence_auto_test, h_convergence_auto_test
export get_refined_models, get_distributed_refined_models, get_octree_refined_models, get_3D_octree_refined_models
export nref, nc, nc_horizontal, dx, dx_horizontal
export convergence_rate

end
