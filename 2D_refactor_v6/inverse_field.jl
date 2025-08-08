struct ForwardMapField{A}  <: Field
  p::A
end


function Gridap.Arrays.return_cache(f::ForwardMapField,
  cellx::AbstractArray{<:VectorValue{2}})

  y = similar(cellx,VectorValue{3,Float64})
  return y
end


function Gridap.Arrays.evaluate!(cache,f::ForwardMapField,
  cellx::AbstractArray{<:VectorValue{2}} )
  y = cache

  p = f.p



  if p == 1
    map!(x -> VectorValue(1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 ),
                          1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 ) * tan(x[1]),
                          1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 ) * tan(x[2]) ),
                                y, cellx)
  elseif p == 2
    map!(x -> VectorValue(1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 )* tan(x[2]),
                          1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 ) * tan(x[1]),
                          1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 )  ),
                                y, cellx)
  elseif p == 3
    map!(x -> VectorValue(1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 )* tan(x[1]),
                          1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 ),
                          1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 ) * tan(x[2])  ),
                                y, cellx)
  elseif p == 4
    map!(x -> VectorValue(-1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 ),
                          -1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 ) * tan(x[1]),
                          1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 ) * tan(x[2]) ),
                                y, cellx)
  elseif p == 5
    map!(x -> VectorValue(-1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 )* tan(x[2]),
                          1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 ) * tan(x[1]),
                          -1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 )  ),
                                y, cellx)
  elseif p == 6
    map!(x -> VectorValue(-1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 )* tan(x[1]),
                          -1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 ),
                          1/sqrt( 1 + (tan(x[1]))^2 + (tan(x[2]))^2 ) * tan(x[2])  ),
                                y, cellx)
  end



  return y
end


function Gridap.Arrays.return_cache(f::ForwardMapField,x::VectorValue{2})
  y = zero(VectorValue{3,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::ForwardMapField,x::VectorValue{2})
  y = cache

  p = f.p

  α,β = x


  rho = sqrt( 1 + (tan(α))^2 + (tan(β))^2  )

  if p == 1
    y[1] = 1/rho
    y[2] = 1/rho * tan(α)
    y[3] = 1/rho * tan(β)
  elseif p == 2
    y[3] = 1/rho
    y[2] = 1/rho * tan(α)
    y[1] = 1/rho * tan(β)
  elseif p == 3
    y[2] = 1/rho
    y[1] = 1/rho * tan(α)
    y[3] = 1/rho * tan(β)
  elseif p == 4
    y[1] = -1/rho
    y[3] = -1/rho * tan(α)
    y[3] = 1/rho * tan(β)
  elseif p == 5
    y[3] = -1/rho
    y[2] = 1/rho * tan(α)
    y[1] = -1/rho * tan(β)
  elseif p == 6
    y[2] = -1/rho
    y[1] = -1/rho * tan(α)
    y[3] = 1/rho * tan(β)
  end


  return y
end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:ForwardMapField},
  cellx::AbstractArray{<:VectorValue{2,T}}) where {T}
  _T = typeof(TensorValue{2,3,T}) ## Jacobian is 3x2, recall we need to return transpose
  y = similar(cellx,_T,size(cellx))
  CachedArray(y)
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:ForwardMapField},cellx::AbstractArray{<:VectorValue{2}})
  cache,  = c
  setsize!(cache,size(cellx))
  y = cache.array
  p = f.object.p

  @check p == 1

  map!(x -> TensorValue{2,3}( -tan(x[1])*(sec(x[1]))^2/( ( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(3/2) ),
                              -tan(x[2])*(sec(x[2]))^2/( ( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(3/2) ),
                            ( -tan(x[1])*(sec(x[1]))^2/( ( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(3/2) ) * tan(x[1])
                                + 1/( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(1/2) * (sec(x[1]))^2   ),
                              -tan(x[2])*(sec(x[2]))^2/( ( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(3/2) ) * tan(x[1]),
                              -tan(x[1])*(sec(x[1]))^2/( ( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(3/2) ) * tan(x[2]),
                            ( -tan(x[2])*(sec(x[2]))^2/( ( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(3/2) ) * tan(x[2])
                            + 1/( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(1/2) * (sec(x[2]))^2 )                   ,
                            ),
                    y, cellx)

  return y

end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:ForwardMapField},x::VectorValue{2,T}) where {T}
  zero(TensorValue{2,3,T}) ## Jacobian is 3x2, recall we need to return transpose
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:ForwardMapField},x::VectorValue{2})
  y = cache

  p = f.object.p

  @check p == 1


  y = TensorValue{2,3}( -tan(x[1])*(sec(x[1]))^2/( ( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(3/2) ),
                        -tan(x[2])*(sec(x[2]))^2/( ( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(3/2) ),
                      ( -tan(x[1])*(sec(x[1]))^2/( ( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(3/2) ) * tan(x[1])
                          + 1/( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(1/2) * (sec(x[1]))^2   ),
                        -tan(x[2])*(sec(x[2]))^2/( ( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(3/2) ) * tan(x[1]),
                        -tan(x[1])*(sec(x[1]))^2/( ( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(3/2) ) * tan(x[2]),
                      ( -tan(x[2])*(sec(x[2]))^2/( ( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(3/2) ) * tan(x[2])
                       + 1/( 1 + (tan(x[1]))^2 + (tan(x[2]))^2  )^(1/2) * (sec(x[2]))^2 )                   ,
                          )
  return y
end



################################################################################
struct InverseMapField{A}  <: Field
  p::A
end


function Gridap.Arrays.return_cache(f::InverseMapField,
  cellx::AbstractArray{<:VectorValue{3}})

  y = similar(cellx,VectorValue{2,Float64})
  return y
end


function Gridap.Arrays.evaluate!(cache,f::InverseMapField,
  cellx::AbstractArray{<:VectorValue{3}} )
  y = cache

  p = f.p

  if p == 1
    map!(x -> VectorValue(atan(x[2],x[1]), atan(x[3],x[1]) ), y, cellx)
  elseif p == 2
    map!(x -> VectorValue(atan(x[2],x[3]), atan(x[1],x[3]) ), y, cellx)
  elseif p == 3
    map!(x -> VectorValue(atan(x[1],x[2]), atan(x[3],x[2]) ), y, cellx)
  elseif p == 4
    map!(x -> VectorValue(atan(-x[2],-x[1]), atan(x[3],-x[1]) ), y, cellx)
  elseif p == 5
    map!(x -> VectorValue(atan(x[2],-x[3]), atan(-x[1],-x[3]) ), y, cellx)
  elseif p == 6
    map!(x -> VectorValue(atan(-x[1],-x[2]), atan(x[3],-x[2]) ), y, cellx)
  end

  return y
end


function Gridap.Arrays.return_cache(f::InverseMapField,x::VectorValue{3})
  y = zero(VectorValue{2,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::InverseMapField,x::VectorValue{3})
  X,Y,Z = x

  p = f.p

  y = cache

  if p == 1
    y[1] = atan(Y,X)
    y[2] = atan(Z,X)
  elseif p == 2
    y[1] = atan(Y,Z)
    y[2] = atan(X,Z)
  elseif p == 3
    y[1] = atan(X,Y)
    y[2] = atan(Z,Y)
  elseif p == 4
    y[1] = atan(-Y,-X)
    y[2] = atan(Z,-X)
  elseif p == 5
    y[1] = atan(Y,-Z)
    y[2] = atan(-X,-Z)
  elseif p == 6
    y[1] = atan(-X,-Y)
    y[2] = atan(Z,-Y)
  end


  return y
end
