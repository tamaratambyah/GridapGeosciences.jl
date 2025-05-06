"""
MetricMeasure

Applies the measure to an array of metrics
"""
struct MetricMeasure <: Map end

function Gridap.Arrays.return_cache(k::MetricMeasure,gx::AbstractVector{<:TensorValue})
  g = first(gx)
  T = typeof( sqrt(meas(g)) )
  y = similar(gx,T)
  return y
end

function Gridap.Arrays.evaluate!(cache,k::MetricMeasure,gx::AbstractVector{<:TensorValue})
  y = cache
  map!(x -> sqrt(meas(x)), y, gx)
  return y
end

"""
MetricInverse

Applies the inverse to an array of metrics
"""
struct MetricInverse <: Map end

function Gridap.Arrays.return_cache(k::MetricInverse,gx::AbstractVector{<:TensorValue})
  g = first(gx)
  T = typeof( inv(g) )
  y = similar(gx,T)
  return y
end

function Gridap.Arrays.evaluate!(cache,k::MetricInverse,gx::AbstractVector{<:TensorValue})
  y = cache
  map!(x -> inv(x), y, gx)
  return y
end


"""
LazyMult

Multiplies two lazy arrays
"""
struct LazyMult <: Map end

function Gridap.Arrays.return_cache(k::LazyMult,ax::AbstractArray,bx::AbstractArray)
  @check length(ax) == length(bx)
  a = first(ax)
  b = first(bx)
  T = typeof( a*b )
  y = similar(bx,T)
  return y
end

function Gridap.Arrays.evaluate!(cache,k::LazyMult,ax::AbstractArray,bx::AbstractArray)
  y = cache
  y .= ax .* bx
  return y
end
