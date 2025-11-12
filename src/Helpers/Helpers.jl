module  Helpers

using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using Gridap.Algebra, Gridap.FESpaces, Gridap.Helpers
import Gridap.TensorValues: MultiValue
using LinearAlgebra
using FillArrays


include("forward_map_2D.jl")
include("forward_map_3D.jl")
include("inverse_map.jl")
include("ForwardMap.jl")
include("overloads.jl")
include("analytical_functions.jl")
include("analytical_functions_autodiff.jl")
include("operators.jl")
include("coordinate_mappings.jl")
include("vector_projection_analytic_functions.jl")


export ForwardMap
export forward_map_2D, forward_jacobian_2D, covarient_basis, forward_pinv_jacobian
export forward_map_3D, forward_jacobian_3D, forward_jacobian
export contravariant_basis

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
export g_star
export xyz2θϕr, θϕ2xyz, spherical_to_cartesian_matrix, xyz2θϕ

export forward_maps, RADIUS, THICKNESS

end
