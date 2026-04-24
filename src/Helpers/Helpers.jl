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
# whole GridapGeoosciences package in this Helpers module, without creating circular dependencies
using GridapGeosciences.Fields
import GridapGeosciences.Fields: J

include("overloads.jl")
include("operators.jl")
include("coordinate_mappings.jl")
include("vector_projection_analytic_functions.jl")

export forward_jacobian, covariant_basis

export pinvJ
export sqrtg,  detg
export grad_meas
export metric, inv_metric
export surflap, surfdiv, contr_gradf, sgrad

export panel_to_cartesian
export normal_vec, tangent_vec
export contra_v, projection_v
export piola
export xyz2θϕr, θϕ2xyz, xyz2θϕ


end
