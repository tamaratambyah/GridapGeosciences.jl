module  Helpers

using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using Gridap.Algebra, Gridap.FESpaces
using LinearAlgebra
using FillArrays

import ..RADIUS

include("forward_map.jl")
include("inverse_map.jl")
# include("Helpers/overloads.jl")
include("analytical_functions.jl")
include("coordinate_mappings.jl")
include("vector_projection_analytic_functions.jl")

export forward_map, forward_jacobian, covarient_basis, forward_pinv_jacobian
export contravariant_basis

# export f_sin,fθϕ,f_XYZ,fX
export rho, rho3
export sqrtg, _sqrtg, detg
export grad_meas
export analytic_metric, analytic_inv_metric, _analytic_inv_metric
# export analytic_J1
export analytic_perp_matrix
export surflap, surfdiv, contr_gradf, sgrad

export panel_to_cartesian, panel_to_latlon
export normal_vector_from_basis
export vector_length, normal_vec, tangent_vec
export contra_v, contra_v_comp, contra_v_perp, projection_v

export xyz2θϕr, θϕ2xyz, spherical_to_cartesian_matrix, xyz2θϕ

export E, F, G

end
