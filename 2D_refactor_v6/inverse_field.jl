using Gridap.Helpers
struct InverseMapField2{A}  <: Field
  p::A
end


function Gridap.Arrays.return_cache(f::InverseMapField2,
  cellx::AbstractArray{<:VectorValue{3}})

  y = similar(cellx,VectorValue{2,Float64})
  return y
end


function Gridap.Arrays.evaluate!(cache,f::InverseMapField2,
  cellx::AbstractArray{<:VectorValue{3}} )
  y = cache

  p = f.p

  if p == 1
    map!(x -> VectorValue(atan(x[2],x[1]), atan(x[3],x[1]) ), y, cellx)
  elseif p == 6
    map!(x -> VectorValue(atan(-x[3],x[2]), atan(x[1],x[2]) ), y, cellx)
  end

  return y
end


function Gridap.Arrays.return_cache(f::InverseMapField2,x::VectorValue{3})
  y = zero(VectorValue{2,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::InverseMapField2,x::VectorValue{3})
  X,Y,Z = x

  p = f.p

  y = cache
  if p == 1
    y[1] = atan(Y,X)
    y[2] = atan(Z,X)
  elseif p == 6
    y[1] = atan(-Z,Y)
    y[2] = atan(X,Y)
  end


  return y
end
