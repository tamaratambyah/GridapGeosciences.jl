"""
Bump
- BumpMap
- BumpField

Map points on the reference panel (panel 1) from 2D <-> 3D.

The `A_bump`, `B_bump`, `b_bump` Tensor(Vector)Values are based on the rotation
of panel 1 in the C-CAM ordering

The dimension of cellx is used to determine what map to apply via type dispatching:
  if D == 3, bump 3D -> 2D
  if D == 2, bump 2D -> 3D

The same functionality is provided for a GridapMap and GridapField
- the GridapField is used to map cell coordinates and cell maps
    * the gradient is implemented
- the GridapMap is generally not used, but is provided regardless for testing
purposes

This map is the inverse of itself
"""

"""
BumpField
"""

"""
BumpMap
"""
struct BumpField{A,B,b} <: Field
  A_bump::A
  B_bump::B
  b_bump::b
end

"""
D == 3, -> y == 2 components; bump 3D -> 2D
y .= A.⋅cellx
"""
function Gridap.Arrays.return_cache(f::BumpField,
  cellx::AbstractArray{<:VectorValue{3}})
  A = f.A_bump
  x = first(cellx)

  T = typeof(A⋅x)
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,f::BumpField,
  cellx::AbstractArray{<:VectorValue{3}})
  y = cache
  A = f.A_bump
  map!(x -> A⋅x, y, cellx)
  return y
end


function Gridap.Arrays.return_cache(f::BumpField,x::VectorValue{3})
  A = f.A_bump
  T = typeof(A⋅x)
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::BumpField,x::VectorValue{3})
  y = cache
  A = f.A_bump
  y = A.⋅x
  return y
end

### gradient is just the coefficient matrix
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:BumpField},
  cellx::AbstractArray{<:VectorValue{3}})
  T = typeof(transpose(f.object.A_bump) )
  y = similar(cellx,T)
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:BumpField},
  cellx::AbstractArray{<:VectorValue{3}})
  setsize!(cache,size(cellx))
  y = cache.array
  AT = transpose(f.object.A_bump)
  map!(x -> AT, y, cellx)

  return y
end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:BumpField},x::VectorValue{3})
  T = typeof(transpose(f.object.A_bump) )
  y = zero(T)

  return y
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:BumpField},x::VectorValue{3})
  y = cache
  y = transpose(f.object.A_bump)

  return y
end



"""
D == 2 -> y == 3 components;  bump 2D -> 3D
y .= B.⋅cellx .+ b
"""
function Gridap.Arrays.return_cache(f::BumpField,
  cellx::AbstractArray{<:VectorValue{2}})
  B = f.B_bump
  x = first(cellx)

  T = typeof(B⋅x)
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,f::BumpField,
  cellx::AbstractArray{<:VectorValue{2}})
  y = cache
  B = f.B_bump
  b = f.b_bump

  map!(x -> B⋅x+b, y, cellx)
  return y
end


function Gridap.Arrays.return_cache(f::BumpField,x::VectorValue{2})
  B = f.B_bump
  T = typeof(B⋅x)
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::BumpField,x::VectorValue{2})
  y = cache
  B = f.B_bump
  b = f.b_bump

  y = B.⋅x .+ b
  return y
end

### gradient is just the coefficient matrix
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:BumpField},
  cellx::AbstractArray{<:VectorValue{2}})
  T = typeof(transpose(f.object.B_bump) )
  y = similar(cellx,T)
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:BumpField},
  cellx::AbstractArray{<:VectorValue{2}})
  setsize!(cache,size(cellx))
  y = cache.array
  AT = transpose(f.object.B_bump)
  map!(x -> AT, y, cellx)

  return y
end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:BumpField},x::VectorValue{2})
  T = typeof(transpose(f.object.B_bump) )
  y = zero(T)

  return y
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:BumpField},x::VectorValue{2})
  y = cache
  y = transpose(f.object.B_bump)

  return y
end



"""
BumpMap
"""
struct BumpMap{A,B,b} <: Map
  A_bump::A
  B_bump::B
  b_bump::b
end

BumpMap() = BumpMap(A_bump,B_bump,b_bump)


"""
D == 3, -> y == 2 components; bump 3D -> 2D
y .= A.⋅cellx
"""
function Gridap.Arrays.return_cache(f::BumpMap,
  cellx::AbstractArray{<:VectorValue{3}})
  A = f.A_bump
  x = first(cellx)

  T = typeof(A⋅x)
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,f::BumpMap,
  cellx::AbstractArray{<:VectorValue{3}})
  y = cache
  A = f.A_bump
  map!(x -> A⋅x, y, cellx)
  return y
end


function Gridap.Arrays.return_cache(f::BumpMap,x::VectorValue{3})
  A = f.A_bump
  T = typeof(A⋅x)
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::BumpMap,x::VectorValue{3})
  y = cache
  A = f.A_bump
  y = A.⋅x
  return y
end

"""
D == 2 -> y == 3 components;  bump 2D -> 3D
y .= B.⋅cellx .+ b
"""
function Gridap.Arrays.return_cache(f::BumpMap,
  cellx::AbstractArray{<:VectorValue{2}})
  B = f.B_bump
  x = first(cellx)

  T = typeof(B⋅x)
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,f::BumpMap,
  cellx::AbstractArray{<:VectorValue{2}})
  y = cache
  B = f.B_bump
  b = f.b_bump

  map!(x -> B⋅x+b, y, cellx)
  return y
end


function Gridap.Arrays.return_cache(f::BumpMap,x::VectorValue{2})
  B = f.B_bump
  T = typeof(B⋅x)
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::BumpMap,x::VectorValue{2})
  y = cache
  B = f.B_bump
  b = f.b_bump

  y = B.⋅x .+ b
  return y
end



"""
bump_matrices

returns the TensorValue matrices required to bump points 2D <-> 3D on panel 1.
note, it is assumed the cube faces are [-a,a]^2 so that panel 1 has X = a
"""

function bump_matrics(a::Float64)
  _A , _B, _b = _bump_matrics(a)
  A_bump = TensorValue(_A)
  B_bump = TensorValue(_B)
  b_bump = VectorValue(_b)
  return A_bump, B_bump, b_bump
end

function _bump_matrics(a::Float64)
  A = [0.0 1.0 0.0
       0.0 0.0 1.0]
  B = [0.0 0.0
       1.0 0.0
       0.0 1.0]
  b = a .* [1.0
                0.0
                0.0]
  return A, B, b
end

const A_bump, B_bump, b_bump = bump_matrics(a)
