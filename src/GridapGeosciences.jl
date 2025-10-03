module GridapGeosciences
using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers, Gridap.Visualization
using Gridap.Algebra, Gridap.FESpaces
using LinearAlgebra
using FillArrays

const global RADIUS = 1.0
export RADIUS

include("Helpers/Helpers.jl")

include("Fields/Fields.jl")

include("Geometry/Geometry.jl")

include("Adaptivity/Adaptivity.jl")

include("ODEs/ODEs.jl")

include("Visualisation/Visualisation.jl")

include("Distributed/Distributed.jl")

include("Exports.jl")

end
