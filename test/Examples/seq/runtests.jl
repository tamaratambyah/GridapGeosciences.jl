using Test

@testset "AdvectionTutorial" begin include("AdvectionTest.jl") end
@testset "LaplaceBeltramiTutorial" begin include("../LaplaceBeltrami.jl") end
