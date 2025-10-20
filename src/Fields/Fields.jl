module Fields

using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using Gridap.Algebra, Gridap.FESpaces
using LinearAlgebra
using FillArrays

import ..RADIUS

using GridapGeosciences.Helpers
import GridapGeosciences.Helpers: dXda, dXdb, dYda,dYdb, dZda, dZdb
import GridapGeosciences.Helpers: forward_map, inverse_map, inverse_jacobian

include("ForwardMap.jl")
include("InverseMap.jl")
include("MatMultField.jl")
include("AffineField.jl")
include("Cartesian2SphericalMap.jl")

export ForwardMapPanel1
export MatMultField, MyAffineField
export Cartesian2SphereicalMap


end
