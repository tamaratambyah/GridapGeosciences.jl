using GridapGeosciences
using Test

TESTCASE = get(ENV, "TESTCASE", "seq")

# Sequential tests
if TESTCASE ∈ ("all", "seq", "seq-l2-projection")
  include("Projection/seq/runtests.jl")
end

if TESTCASE ∈ ("all", "seq", "seq-laplacian")
  include("Laplacian/seq/runtests.jl")
end


# MPI tests
if TESTCASE ∈ ("all", "mpi", "mpi-l2-projection")
   include("Projection/mpi/runtests.jl")
end

if TESTCASE ∈ ("all", "mpi", "mpi-laplacian")
  include("Laplacian/mpi/runtests.jl")
end
