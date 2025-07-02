
struct ShiftField{A}  <: Field
  shifts::A
end

"""
D == 3, -> y == 2 components; bump 3D -> 2D
y .= A.⋅cellx
"""
function Gridap.Arrays.return_cache(f::ShiftField,
  cellx::AbstractArray{<:VectorValue{2}})

  x = first(cellx)
  T = typeof(x)
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,f::ShiftField,
  cellx::AbstractArray{<:VectorValue{2}} )
  y = cache
  v = f.shifts
  map!(x -> x + v, y, cellx)


  return y
end


function Gridap.Arrays.return_cache(f::ShiftField,x::VectorValue{2})
  T = typeof(x)
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::ShiftField,x::VectorValue{2})
  v = f.shifts
  y = cache

  y = x + v

  return y
end

### gradient is just the coefficient matrix
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:ShiftField},
  cellx::AbstractArray{<:VectorValue{2}})
  _T = typeof(TensorValue{2,2,Float64})
  y = similar(cellx,_T,size(cellx))
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:ShiftField},
  cellx::AbstractArray{<:VectorValue{2}})
  cache, = c
  setsize!(cache,size(cellx))
  y = cache.array

  map!(x -> one((TensorValue{2,2,Float64})), y, cellx)

  return y
end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:ShiftField},x::VectorValue{2})
  zero(TensorValue{2,2,Float64})

end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:ShiftField},x::VectorValue{2})
  y = cache
  y = one(TensorValue{2,2,Float64})

  return y
end


###############################################################################

struct InversionField{A}  <: Field
  P::A
end

"""
D == 3, -> y == 2 components; bump 3D -> 2D
y .= A.⋅cellx
"""
function Gridap.Arrays.return_cache(f::InversionField,
  cellx::AbstractArray{<:VectorValue{2}})

  x = first(cellx)
  T = typeof(x)
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,f::InversionField,
  cellx::AbstractArray{<:VectorValue{2}} )
  y = cache
  P =  f.P
  map!(x -> P ⋅ x, y, cellx)

  return y
end


function Gridap.Arrays.return_cache(f::InversionField,x::VectorValue{2})
  T = typeof(x)
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::InversionField,x::VectorValue{2})
  P =  f.P
  y = cache

   y = P .⋅ x

  return y
end

### gradient is just the coefficient matrix
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:InversionField},
  cellx::AbstractArray{<:VectorValue{2}})
  _T = typeof(TensorValue{2,2,Float64})
  y = similar(cellx,_T,size(cellx))
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:InversionField},
  cellx::AbstractArray{<:VectorValue{2}})
  setsize!(cache,size(cellx))
  y = cache.array
  P =  f.object.P
  map!(x -> transpose(P), y, cellx)

  return y
end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:InversionField},x::VectorValue{2})
  zero(TensorValue{2,2,Float64})

end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:InversionField},x::VectorValue{2})
  y = cache
  P =  f.object.P
  y = transpose(P)

  return y
end
