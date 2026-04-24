module Fields

using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using Gridap.Algebra, Gridap.FESpaces
using LinearAlgebra
using FillArrays


using GridapGeosciences.Helpers
import GridapGeosciences.Helpers: forward_jacobian

include("MatMultField.jl")
include("AffineField.jl")
include("Cartesian2SphericalMap.jl")
include("Cartesian2SphericalMap3D.jl")

export MatMultField, MyAffineField
export Cartesian2SphericalMap
export Cartesian2SphericalMap3D

end
