using Test

@testset "AmbientHodgeLaplacianScalarTests" begin include("AmbientHodgeLaplacianScalarTests.jl") end

@testset "AmbientLaplaceBeltrami" begin include("AmbientLaplaceBeltramiTests.jl") end

@testset "AmbientLinearisedShallowWater" begin include("AmbientLinearisedShallowWaterTests.jl") end

@testset "AmbientSurfaceArea" begin include("AmbientSurfaceAreaTests.jl") end
