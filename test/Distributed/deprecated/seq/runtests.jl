using Test

@testset "DistributedNormalVectorTests" begin include("DistributedNormalVectorTests.jl") end

@testset "DistributedNormalVectorTests3D" begin include("DistributedNormalVectorTests3D.jl") end

@testset "DistributedPanelIds" begin include("DistributedPanelIdsTests.jl") end

@testset "DistributedSurfaceArea" begin include("DistributedSurfaceAreaTests.jl") end

@testset "DistributedTriangulationTests" begin include("DistributedTriangulationTests.jl") end
