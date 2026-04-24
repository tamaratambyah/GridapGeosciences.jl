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


using GridapGeosciences.Fields
import GridapGeosciences.Fields: ForwardMap2D, ForwardMap3D, ForwardMap2Dor3D

include("overloads.jl")
include("analytical_functions_autodiff.jl")
include("operators.jl")
include("coordinate_mappings.jl")
include("vector_projection_analytic_functions.jl")
include("convergence_tools.jl")

export forward_jacobian, covariant_basis

export pinvJ
export sqrtg,  detg
export grad_meas
export metric, inv_metric
export surflap, surfdiv, contr_gradf, sgrad

export panel_to_cartesian, panel_to_latlon
export normal_vec, tangent_vec
export contra_v, projection_v
export piola
export xyz2θϕr, θϕ2xyz, spherical_to_cartesian_matrix, xyz2θϕ

export p_convergence_auto_test, h_convergence_auto_test
export get_refined_models, get_distributed_refined_models, get_octree_refined_models, get_3D_octree_refined_models
export nref, nc, nc_horizontal, dx, dx_horizontal
export convergence_rate

end
