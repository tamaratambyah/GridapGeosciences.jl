module Fields

using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using Gridap.Algebra, Gridap.FESpaces
using LinearAlgebra
using FillArrays


include("ForwardMap.jl")
include("MatMultField.jl")
include("AffineField.jl")
include("Cartesian2SphericalMap.jl")

export ForwardMap, ForwardMap2D, ForwardMap3D, ForwardMap2DGenerator, ForwardMap3DGenerator
export ForwardMap2Dor3D, ForwardMap2Dor3DGenerator
export J
export MatMultField, MyAffineField
export Cartesian2SphericalMap

end
