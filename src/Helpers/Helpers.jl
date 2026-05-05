module  Helpers

using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using Gridap.Algebra, Gridap.FESpaces, Gridap.Helpers, Gridap.Arrays
import Gridap.TensorValues: MultiValue
using LinearAlgebra
using FillArrays

using GridapGeosciences.Fields
import GridapGeosciences.Fields: J, normal_vec

include("Overloads.jl")
include("Operators.jl")
include("CoordinateMappings.jl")
include("VectorPullback.jl")

export forward_jacobian, covariant_basis

export pinvJ
export sqrtg,  detg
export grad_meas
export metric, inv_metric
export surflap, surfdiv, contr_gradf, sgrad

export panel_to_cartesian
export tangent_vec
export contra_v
export piola
export xyz2θϕr, θϕ2xyz, xyz2θϕ


end
