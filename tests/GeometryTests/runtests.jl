module GeometryTests

using Test

@testset "CellMap" begin include("CellMapTests.jl") end

@testset "NormalVector" begin include("NormalVectorTests.jl") end

@testset "Refinement" begin include("RefinementTests.jl") end

@testset "PanelIds" begin include("PanelIdsTests.jl") end

end

# @testset "CellMap" begin include("tests/GeometryTests/CellMapTests.jl") end
