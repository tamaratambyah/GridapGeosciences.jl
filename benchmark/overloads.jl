## We have to overwrite this function for CountFlops to be compatible with
## UnstructuredDiscreteModel. This is because it does not know that evaluate on
## the iterator is type stable (and for some reason does a ccal).
## In this function, provide the type to vals explicitily.
## The length of this compressed array is one, so it should skip the for loop
# function Gridap.Arrays._lazy_map_compressed(g::CompressedArray...)
#   #vals = map(evaluate, map(gi->gi.values,g)...)
#   n = length(first(g).values)
#   v1 = evaluate(map(gi -> first(gi.values),g)...)
#   vals = Vector{typeof(v1)}(undef,n)
#   vals[1] = v1
#   for i in 2:n
#     vals[i] = evaluate(map(gi -> gi.values[i],g)...)
#   end
#   ptrs = first(g).ptrs
#   Gridap.Arrays.CompressedArray(vals,ptrs)
# end
function Gridap.Arrays._lazy_map_compressed(g::CompressedArray...)
  #vals = map(evaluate, map(gi->gi.values,g)...)
  n = length(first(g).values)
  v1 = evaluate(map(gi -> first(gi.values),g)...)
  vals = Vector{typeof(v1)}(undef,n)
  for i in 1:n
    vals[i] = evaluate(map(gi -> gi.values[i],g)...)
  end
  ptrs = first(g).ptrs
  Gridap.Arrays.CompressedArray(vals,ptrs)
end


## We have to overwrite this function for CountFlops to be compatible with
## UnstructuredDiscreteModel. This is because it does not know that evaluate on
## the iterator is type stable (and for some reason does a ccal).
## In this case, N = 1
function Gridap.Arrays.CachedArray(T,N)
  # s = tuple([0 for i in 1:N]...)
  s = tuple([0]...)
  a = Array{T,N}(undef,s)
  Gridap.Arrays.CachedArray(a)
end
