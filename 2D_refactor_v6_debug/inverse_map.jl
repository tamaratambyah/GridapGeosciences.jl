function inverse_map(XYZ,p)
  X,Y,Z = XYZ

  if p == 1
    α = atan(Y,X)
    β = atan(Z,X)
  elseif p == 2
    α = atan(Y,Z)
    β = -atan(X,Z)
  elseif p == 3
    α = -atan(X,Y)
    β = atan(Z,Y)
  end
  Point(α,β)

end

function inverse_jacobian(p)
  function _inverse_jacobian(αβ)
    XYZ = forward_map(αβ,p)
    transpose(inverse_jacobian(XYZ,p))

  end
end

function inverse_jacobian(XYZ,p)
  X,Y,Z = XYZ

  if p == 1
    dadX = - Y/(X^2 + Y^2)
    dadY = X/(X^2 + Y^2)
    dadZ = 0.0
    dbdX = -Z/(X^2 + Z^2)
    dbdY = 0.0
    dbdZ = X/(X^2 + Z^2)
  elseif p == 2
    dadX = 0.0
    dadY = Z/(Y^2 + Z^2)
    dadZ = -Y/(Y^2 + Z^2)
    dbdX = -Z/(X^2 + Z^2)
    dbdY = 0.0
    dbdZ = X/(X^2 + Z^2)
  elseif p == 3
    dadX = -Y/(X^2 + Y^2)
    dadY = X/(X^2 + Y^2)
    dadZ = 0.0
    dbdX = 0.0
    dbdY = -Z/(Y^2 + Z^2)
    dbdZ = Y/(Y^2 + Z^2)
  end

  ## J = [dadX dadY dadX
  ##      dbdX dbdY dbdZ ]
  ## As a TensorValue data = (dadX,dadY,dadZ, dbdX, dbdY, dbdZ)

  TensorValue{2,3}(dadX,dbdX, dadY,dbdY, dadZ,dbdZ)

end



################################################################################
################################################################################
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
