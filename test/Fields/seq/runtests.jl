using Test

@testset "AffineField" begin include("AffineFieldTests.jl") end

@testset "ForwardInverseMap" begin include("ForwardInverseMapTests.jl") end

@testset "Overloads" begin include("OverloadTests.jl") end
