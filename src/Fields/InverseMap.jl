struct InverseMap{A}  <: Field
  p::A
end


function Gridap.Arrays.return_cache(f::InverseMap,cellx::AbstractArray{<:VectorValue{3}})
  p = f.p
  x = first(cellx)

  T = typeof(inverse_map(x,p))
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,f::InverseMap,cellx::AbstractArray{<:VectorValue{3}} )
  y = cache
  p = f.p

  map!(x -> inverse_map(x,p), y, cellx)

  return y
end


function Gridap.Arrays.return_cache(f::InverseMap,x::VectorValue{3})
  p = f.p

  T = typeof(inverse_map(x,p))
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::InverseMap,x::VectorValue{3})
  y = cache

  p = f.p
  y = inverse_map(x,p)

  return y
end

### gradient is just the coefficient matrix
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:InverseMap},
  cellx::AbstractArray{<:VectorValue{3}})

  p = f.object.p
  x = first(cellx)
  T = typeof(transpose(inverse_jacobian(x,p)))

  y = similar(cellx,T)
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:InverseMap},
  cellx::AbstractArray{<:VectorValue{3}})
  cache, = c
  setsize!(cache,size(cellx))

  y = cache
  p = f.object.p
  map!(x -> transpose(inverse_jacobian(x,p)), y, cellx)
  return y

end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:InverseMap},x::VectorValue{3})
  transpose( inverse_jacobian(x,p) )
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:InverseMap},x::VectorValue{3})
  y = cache
  p = f.object.p
  y = transpose(inverse_jacobian(x,p))
  return y
end
