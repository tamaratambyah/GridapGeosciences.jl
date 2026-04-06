module TransientWaveEquationTestsSeq
using PartitionedArrays
include("../TransientWaveEquation.jl")

with_debug() do distribute
  TransientWaveEquation.main_transient(distribute,1;n_ref_lvls=3,return_vtk=true)
end

end

# @testset "TransientWave" begin include("test/Transient/seq/TransientWaveEquationTests.jl") end
