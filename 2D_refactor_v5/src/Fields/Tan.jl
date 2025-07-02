struct TanField{A} <: Field
  a::A
end

function Gridap.Arrays.return_cache(f::TanField,
  cellx::AbstractArray{<:VectorValue{2}})
  x = first(cellx)
  T = typeof(x)
  y = similar(cellx,T)
  return y
end


function Gridap.Arrays.evaluate!(cache,f::TanField,
  cellx::AbstractArray{<:VectorValue{2}})
  y = cache
  a = f.a
  map!(x -> VectorValue(atan(x[1]/a),atan(x[2]/a)), y, cellx)
  return y
end


function Gridap.Arrays.return_cache(f::TanField,x::VectorValue{2})
  T = typeof(x)
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::TanField,x::VectorValue{2})
  y = cache
  a = f.a
  y = VectorValue(atan(x[1]/a),atan(x[2]/a))
  return y
end

### gradient is just the coefficient matrix
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:TanField},
  cellx::AbstractArray{<:VectorValue{2}})
  _T = typeof(TensorValue{2,2,Float64})
  y = similar(cellx,_T,size(cellx))
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:TanField},
  cellx::AbstractArray{<:VectorValue{2}})
  setsize!(cache,size(cellx))
  a = f.object.a
  map!(x -> TensorValue{2,2}( 1/(1+(x[1]/a)^2), 0.0, 0.0, 1/(1+(x[2]/a)^2)
        ),
      y, cellx)

  return y
end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:TanField},x::VectorValue{2})
  zero(TensorValue{2,2,T})

  return y
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:TanField},x::VectorValue{2})
  y = cache
  a = f.object.a
  y = TensorValue{2,2}( 1/(1+(x[1]/a)^2), 0.0, 0.0, 1/(1+(x[2]/a)^2))

  return y
end
