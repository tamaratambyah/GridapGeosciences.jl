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
include("Cartesian2SphereicalMap3D.jl")

export MatMultField, MyAffineField
export Cartesian2SphereicalMap
export Cartesian2SphereicalMap3D

end
