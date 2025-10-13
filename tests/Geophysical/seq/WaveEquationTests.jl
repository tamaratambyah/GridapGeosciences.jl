module WaveEquationTestsSeq
using PartitionedArrays
include("../WaveEquation.jl")

with_debug() do distribute
  WaveEquation.main(distribute,1)
end

end
