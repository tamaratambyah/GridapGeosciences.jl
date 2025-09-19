struct ForwardMapPanel1  <: Field
end

"""
The forward map goes 2D -> 3D.
  X = RADIUS/sqrt(1 + (tan(α))^2 + (tan(β))^2  )
  Y = RADIUS/sqrt(1 + (tan(α))^2 + (tan(β))^2  ) * tan(α)
  Z = RADIUS/sqrt(1 + (tan(α))^2 + (tan(β))^2  ) * tan(β)
"""

function Gridap.Arrays.return_cache(f::ForwardMapPanel1,cellx::AbstractArray{<:VectorValue{2}})
  y = similar(cellx,VectorValue{3,Float64})
  return y
end


function Gridap.Arrays.evaluate!(cache,f::ForwardMapPanel1,cellx::AbstractArray{<:VectorValue{2}} )

  y = cache

  map!(x-> Point(RADIUS/(sqrt(1 + (tan(x[1]))^2 + (tan(x[2]))^2 )),
                 RADIUS/(sqrt(1 + (tan(x[1]))^2 + (tan(x[2]))^2 )) * tan(x[1]),
                 RADIUS/(sqrt(1 + (tan(x[1]))^2 + (tan(x[2]))^2 )) * tan(x[2])  ),
      y, cellx  )

  return y
end


function Gridap.Arrays.return_cache(f::ForwardMapPanel1,x::VectorValue{2})
  y = zero(VectorValue{3,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::ForwardMapPanel1,x::VectorValue{2})
  y = cache
  y = Point(RADIUS/(sqrt(1 + (tan(x[1]))^2 + (tan(x[2]))^2 )),
            RADIUS/(sqrt(1 + (tan(x[1]))^2 + (tan(x[2]))^2 )) * tan(x[1]),
            RADIUS/(sqrt(1 + (tan(x[1]))^2 + (tan(x[2]))^2 )) * tan(x[2])  )

  return y
end

"""
The jacobian is 3 x 2
  J = [dXda dXdb
       dYda dYdb
       dZda dZdb  ]
As a TensorValue data = (dXda,dYda,dZda,   dXdb, dYdb, dZdb)

Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)

The transpose is 2 x 3
  JT = [dXda dYda dZda
        dXdb dYdb dZdb  ]
As a TensorValue data = (dXda,dXdb,  dYda,dYdb,  dZda,dZdb)
"""
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:ForwardMapPanel1},
  cellx::AbstractArray{<:VectorValue{2}})
  y = similar(cellx,TensorValue{2,3,Float64})
  c = CachedArray(y)
  return c
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:ForwardMapPanel1},
  cellx::AbstractArray{<:VectorValue{2}})
  cache, = c
  setsize!(cache,size(cellx))

  y = cache
  map!(x-> TensorValue{2,3}( dXda(x),dXdb(x),  dYda(x),dYdb(x),  dZda(x),dZdb(x) ),
      y, cellx  )
  return y

end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:ForwardMapPanel1},x::VectorValue{2})
  zero(TensorValue{2,3,Float64})
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:ForwardMapPanel1},x::VectorValue{2})
  y = cache

  y = TensorValue{2,3}( dXda(x),dXdb(x),  dYda(x),dYdb(x),  dZda(x),dZdb(x) )

  return y
end
