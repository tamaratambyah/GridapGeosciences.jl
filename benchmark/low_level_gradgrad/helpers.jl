### To save flops, we compute all metric term explicitly
function allg(αβ)

  t1 = tan(αβ[1])
  t2 = tan(αβ[2])

  c1 = cos(αβ[1])
  c2 = cos(αβ[2])

  r = 1.0 + t1*t1 + t2*t2
  r4 = r*r

  c = 1.0/( r4*c1*c1 *c2*c2 )

  f = c * ( -1.0*t1*t2 )
  e = c * ( 1.0 + (t1*t1) )
  g = c * ( 1 + (t2*t2) )

  sqrtg = (e*g - f*f)^(1/2)

  1.0/sqrtg*TensorValue(g, -f, -f, e  )
end

# ### make sure all array caches evaluate
# function Gridap.Arrays.array_cache(dict::Dict,a::Gridap.Arrays.LazyArray)
#   Gridap.Arrays._array_cache!(dict,a)
# end

### benchmark using a lazy collect
function lazy_collect(cache,arr)
  s = 0.0
  for i in eachindex(arr)
    ai = getindex!(cache,arr,i)
    s += ai[1,1]
  end
  return s
end



function total_counts(flops::GFlops.Counter)
  total = 0
  for pn in propertynames(flops)
      total += getproperty(flops, pn)
  end
  return total
end
