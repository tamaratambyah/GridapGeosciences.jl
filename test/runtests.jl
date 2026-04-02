module RunTests

using Test

@testset "Geometry" begin include("GeometryTests/runtests.jl") end

end

# @testset "gs" begin include("test/runtests.jl") end
