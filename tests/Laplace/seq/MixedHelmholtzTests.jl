module MixedHelmholtzTestsSeq
using PartitionedArrays
include("../MixedHelmholtz.jl")

with_debug() do distribute
  MixedHelmholtz.main(distribute,1)
end

end
