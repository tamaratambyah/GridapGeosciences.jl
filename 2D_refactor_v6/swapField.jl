using Gridap.Helpers
struct SwapField{A}  <: Field
  shifts::A
end


function Gridap.Arrays.return_cache(f::SwapField,
  cellx::AbstractArray{<:VectorValue{D}}) where {D}

  A = f.shifts
  x = first(cellx)
  T = typeof(A⋅x)
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,f::SwapField,
  cellx::AbstractArray{<:VectorValue{D}} ) where {D}
  y = cache
  A = f.shifts

  @check D == size(A)[2]
  map!(x -> A ⋅ x, y, cellx)


  return y
end


function Gridap.Arrays.return_cache(f::SwapField,x::VectorValue)
  T = typeof(x)
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SwapField,x::VectorValue{D}) where {D}
  A = f.shifts
  y = cache

  @check D == size(A)[2]

  y = A .⋅ x

  return y
end

### gradient is just the coefficient matrix
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:SwapField},
  cellx::AbstractArray{<:VectorValue})
  A = f.object.shifts
  T = typeof(transpose(A))
  y = similar(cellx,T)
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:SwapField},
  cellx::AbstractArray{<:VectorValue})
  cache, = c
  setsize!(cache,size(cellx))
  y = cache
  A = f.object.shifts
  map!(x -> transpose(A), y, cellx)
  return y

end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:SwapField},x::VectorValue)
  transpose( f.object.shifts )
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:SwapField},x::VectorValue)
  w = cache
  w = f.object.shifts
  return transpose(w)
end
