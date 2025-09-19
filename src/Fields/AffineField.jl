
# Ax + b: multiplication by an invertible matrix and translation of a vector
struct MyAffineField{A,B}  <: Field
  matrix::A
  vector::B
end


function Gridap.Arrays.return_cache(f::MyAffineField,cellx::AbstractArray{<:VectorValue{D}}) where {D}

  A = f.matrix
  b = f.vector
  @check D == size(A)[2]

  x = first(cellx)
  T = typeof(A⋅x + b)
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,f::MyAffineField,cellx::AbstractArray{<:VectorValue{D}} ) where {D}
  y = cache
  A = f.matrix
  b = f.vector

  map!(x -> A ⋅ x + b, y, cellx)


  return y
end


function Gridap.Arrays.return_cache(f::MyAffineField,x::VectorValue{D}) where {D}
  A = f.matrix
  b = f.vector
  @check D == size(A)[2]

  T = typeof(A⋅x + b)
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::MyAffineField,x::VectorValue{D}) where {D}
  y = cache

  A = f.matrix
  b = f.vector
  y = A .⋅ x .+ b

  return y
end

### gradient is just the coefficient matrix
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:MyAffineField},
  cellx::AbstractArray{<:VectorValue})
  A = f.object.matrix
  T = typeof(transpose(A))
  y = similar(cellx,T)
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:MyAffineField},
  cellx::AbstractArray{<:VectorValue})
  cache, = c
  setsize!(cache,size(cellx))
  y = cache
  A = f.object.matrix
  map!(x -> transpose(A), y, cellx)
  return y

end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:MyAffineField},x::VectorValue)
  transpose( f.object.matrix )
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:MyAffineField},x::VectorValue)
  w = cache
  w = f.object.matrix
  return transpose(w)
end
