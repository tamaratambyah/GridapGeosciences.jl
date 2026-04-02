using JLD2
using MPI
using PartitionedArrays
using MPIPreferences

# MPI.Init()
# nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
# ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

# f = jldopen("example.jld2", "r")

nprocs = 2
ranks = with_mpi() do distribute
  distribute(LinearIndices((nprocs,)))
end

a = map(ranks) do r
  r*collect(1:6)
end

map(a) do a
  println(a)
end



map(a,ranks) do a, r
  println("$r")
  jldopen("example_$r.jld2", "a+") do f
    f["$r"] = a
  end
end



_a = map(ranks) do r
  _A = jldopen("example_$r.jld2", "r") do _f
      _f["$r"]
    end
  return _A
end

map(_a) do a
  println(a)
end



struct PVectorLocalStorage{A,B}
  values::A
  indices::B
end

row_partition = uniform_partition(ranks,(1,2),(6,6))
a = pones(row_partition)
