function cube_to_αβ(p)
  if p == 1
    A = [0 1 0
        0 0 1]
  elseif p == 2
    A = [0 1 0
        -1 0 0]
  elseif p == 3
    A = [-1 0 0
          0 0 1]
  elseif p == 4
    A = [0 0 1
        0 1 0]
  elseif p == 5
    A = [-1 0 0
        0 1 0]
  elseif p == 6
    A = [0 0 1
        -1 0 0]
  end

  return TensorValue(A)

end

function rotation_mats(p)
  if p == 1
    A = [1 0 0
         0 1 0
         0 0 1]
  elseif p == 2
    A = [0 0 -1
         0 1 0
         1 0 0]
  elseif p == 3
    A = [0 -1 0
         1 0 0
         0 0 1]
  elseif p == 4
    A = [-1 0 0
          0 0 1
          0 1 0]
  end
  TensorValue(A)
end

using Gridap.Helpers
# left multiply vector by matrix
# works in any dimension, as long as size(A)[2] == length(vec)
struct MatMultField{A}  <: Field
  matrix::A
end


function Gridap.Arrays.return_cache(f::MatMultField,cellx::AbstractArray{<:VectorValue{D}}) where {D}

  A = f.matrix
  @check D == size(A)[2]

  x = first(cellx)
  T = typeof(A⋅x)
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,f::MatMultField,cellx::AbstractArray{<:VectorValue{D}} ) where {D}
  y = cache
  A = f.matrix

  @check D == size(A)[2]
  map!(x -> A ⋅ x, y, cellx)


  return y
end


function Gridap.Arrays.return_cache(f::MatMultField,x::VectorValue{D}) where {D}
  A = f.matrix
  @check D == size(A)[2]

  T = typeof(A⋅x)
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::MatMultField,x::VectorValue{D}) where {D}
  y = cache

  A = f.matrix
  @check D == size(A)[2]

  y = A .⋅ x

  return y
end

### gradient is just the coefficient matrix
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:MatMultField},
  cellx::AbstractArray{<:VectorValue})
  A = f.object.matrix
  T = typeof(transpose(A))
  y = similar(cellx,T)
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:MatMultField},
  cellx::AbstractArray{<:VectorValue})
  cache, = c
  setsize!(cache,size(cellx))
  y = cache
  A = f.object.matrix
  map!(x -> transpose(A), y, cellx)
  return y

end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:MatMultField},x::VectorValue)
  transpose( f.object.matrix )
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:MatMultField},x::VectorValue)
  w = cache
  w = f.object.matrix
  return transpose(w)
end
