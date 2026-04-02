using Test

@testset "LaplaceBeltramiTests" begin include("LaplaceBeltramiTests.jl") end
@testset "HelmholtzTests" begin include("HelmholtzTests.jl") end
@testset "MixedHelmholtzTests" begin include("MixedHelmholtzTests.jl") end

# include("test/Laplace/seq/runtests.jl")
