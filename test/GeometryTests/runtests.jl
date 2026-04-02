module GeometryTests

using Test

@testset "CellMap" begin include("CellMapTests.jl") end

@testset "NormalVector" begin include("NormalVectorTests.jl") end

@testset "Refinement" begin include("RefinementTests.jl") end

@testset "PanelIds" begin include("PanelIdsTests.jl") end

@testset "SurfaceArea" begin include("SurfaceAreaTests.jl") end

end

# @testset "CellMap" begin include("test/GeometryTests/CellMapTests.jl") end
