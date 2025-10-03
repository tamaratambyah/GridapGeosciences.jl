module HelmholtzTestsSeq
using PartitionedArrays
include("../Helmholtz.jl")

with_debug() do distribute
  Helmholtz.main(distribute,1)
end

end
