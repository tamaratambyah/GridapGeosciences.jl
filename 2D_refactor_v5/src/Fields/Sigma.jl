"""
Sigma
- SigmaMap
- SigmaField

The standard spherical mapping from lat-lon <-> 3D Cartesian on sphere surface.

The `latlon` variable is latlon = (θ, ϕ).
θ is the longitude, computed to be [0,2*π] (use rem2pi(θ,RoundDown))
ϕ is the latitude, computed to be [-π/2, π/2]

`r` is the radius of the sphere, which is assumed to be constant

This map is the inverse of itself, dispatch based on number of components
  i.e. cellx = VectorValue{3}, latlon = VectorValue{2}

The same functionality is provided for a GridapMap and GridapField
- the GridapField is used to map cell coordinates and cell maps
    * the gradient is implemented for latlon -> 3D Cartesian on sphere
- the GridapMap is generally not used, but is provided regardless for testing
purposes
"""

"""SigmaField
"""

struct InvSigmaField{A} <: Gridap.Fields.Field
  r::A # sphere radius
end


""" map 3D point on sphere -> latlon (2D), σ^{-1}: (X_s, Y_s, Z_s) → (θ,ϕ)  """
function Gridap.Arrays.return_cache(k::InvSigmaField,cellx::AbstractArray{<:VectorValue{3,T}}) where {T}
  y = similar(cellx,VectorValue{2,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::InvSigmaField,cellx::AbstractArray{<:VectorValue{3}})
  # 3D point on sphere -> lat lon
  y = cache
  map!(x -> VectorValue( rem2pi(atan(x[2], x[1]),RoundDown),
                              atan(x[3], sqrt(x[1]*x[1] + x[2]*x[2]) )
  ),
   y, cellx)
  return y
end

function Gridap.Arrays.return_cache(k::InvSigmaField,x::VectorValue{3,T}) where {T}
  y = zero(VectorValue{2,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::InvSigmaField,x::VectorValue{3})
  # 3D point on sphere -> lat lon
  y = cache
  y = VectorValue( rem2pi(atan(x[2], x[1]),RoundDown),
                  atan(x[3], sqrt(x[1]*x[1] + x[2]*x[2]) )
                   )

  return y
end

"""
gradient of σ^{-1}: (X_s, Y_s, Z_s) → (θ,ϕ)  is:
  J = [ -Y/(X^2+Y^2)                  X/(X^2+Y^2)                0
        -(XZ)/(R^2*sqrt(X^2+Y^2))    -(YZ)/(R^2*sqrt(X^2+Y^2))  sqrt(X^2+Y^2)/R^2 ]

# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
# Note  TensorValue{2,3}(0,4,1,0,5,1) == [0 1 5
                                          4 0 1]

The transpose is:
  JT = [ -Y/(X^2+Y^2)  -(XZ)/(R^2*sqrt(X^2+Y^2))
          X/(X^2+Y^2)  -(YZ)/(R^2*sqrt(X^2+Y^2))
          0             sqrt(X^2+Y^2)/R^2  ]

To write as a TensorValue:
  TensorValue{3,2}( -Y/(X^2+Y^2),   X/(X^2+Y^2),  0,  -(XZ)/(R^2*sqrt(X^2+Y^2)), -(YZ)/(R^2*sqrt(X^2+Y^2)), sqrt(X^2+Y^2)/R^2 )
"""
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:InvSigmaField},
  cellx::AbstractArray{<:VectorValue{3,T}}) where {T}
  _T = typeof(TensorValue{3,2,T})
  y = similar(cellx,_T,size(cellx))
  CachedArray(y)
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:InvSigmaField},cellx::AbstractArray{<:VectorValue{3}})
  cache,  = c
  setsize!(cache,size(cellx))
  y = cache.array
  r = f.object.r

  map!(x -> TensorValue{3,2}( (-x[2]/(x[1]*x[1] + x[2]*x[2])),
                              (x[1]/(x[1]*x[1] + x[2]*x[2])),
                              0.0,
                              (-x[1]*x[3])/(r^2*(sqrt(x[1]*x[1] + x[2]*x[2]))),
                              (-x[2]*x[3])/(r^2*(sqrt(x[1]*x[1] + x[2]*x[2]))),
                              (sqrt(x[1]*x[1] + x[2]*x[2]))/(r^2)
                            ),
                    y, cellx)

  return y

end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:InvSigmaField},x::VectorValue{3,T}) where {T}
  zero(TensorValue{3,2,T}) ## Jacobian is 2x3, recall we need to return transpose
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:InvSigmaField},x::VectorValue{3})
  y = cache
  r = f.object.r
  y = TensorValue{3,2}( (-x[2]/(x[1]*x[1] + x[2]*x[2])),
                         (x[1]/(x[1]*x[1] + x[2]*x[2])),
                         0.0,
                         (-x[1]*x[3])/(r^2*(sqrt(x[1]*x[1] + x[2]*x[2]))),
                         (-x[2]*x[3])/(r^2*(sqrt(x[1]*x[1] + x[2]*x[2]))),
                         (sqrt(x[1]*x[1] + x[2]*x[2]))/(r^2)
                          )
  return y
end





struct SigmaField{A} <: Gridap.Fields.Field
  r::A # sphere radius
end


""" map latlon -> 3D point on sphere, σ: (θ,ϕ) → (X_s, Y_s, Z_s) """
function Gridap.Arrays.return_cache(k::SigmaField,latlon::AbstractArray{<:VectorValue{2,T}}) where {T}
  y = similar(latlon,VectorValue{3,T})
  return y#CachedArray(y)
end

function Gridap.Arrays.evaluate!(cache,f::SigmaField,latlon::AbstractArray{<:VectorValue{2}})
  # lat-lon -> 3D point on sphere
  # latlon = θ, ϕ
  # setsize!(cache,size(cellx))
  y = cache
  r = f.r
  map!(x -> VectorValue(r*cos( rem2pi(x[1],RoundDown) )*cos(x[2]),
                        r*sin( rem2pi(x[1],RoundDown) )*cos(x[2]),
                        r*sin(x[2])),
                        y, latlon)
  return y
end

function Gridap.Arrays.return_cache(k::SigmaField,x::VectorValue{2,T}) where {T}
  y = zero(VectorValue{3,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SigmaField,x::VectorValue{2})
  # lat-lon -> 3D point on sphere
  # latlon = θ, ϕ

  y = cache
  r = f.r
  y = VectorValue(r*cos( rem2pi(x[1],RoundDown) )*cos(x[2]),
                  r*sin( rem2pi(x[1],RoundDown) )*cos(x[2]),
                  r*sin(x[2]))
  return y
end

"""
gradient of σ: θ,ϕ → X_s, Y_s, Z_s  is:
  J = ( -r*sinθ*cosϕ  -r*cosθ*sinϕ
         r*cosθ*cosϕ  -r*sinθ*sinϕ
          0            r*cosϕ  )
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
# Note  TensorValue{2,3}(0,4,1,0,5,1) == [0 1 5
                                          4 0 1]

The transpose is:
  JT = ( -r*sinθ*cosϕ  r*cosθ*cosϕ  0
         -r*cosθ*sinϕ -r*sinθ*sinϕ  r*cosϕ  )

To write as a TensorValue:
  TensorValue{2,3}(-r*sinθ*cosϕ, -r*cosθ*sinϕ, r*cosθ*cosϕ, -r*sinθ*sinϕ, 0, r*cosϕ )
"""
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:SigmaField},
  cellx::AbstractArray{<:VectorValue{2,T}}) where {T}
  _T = typeof(TensorValue{2,3,T})
  y = similar(cellx,_T,size(cellx))
  CachedArray(y)
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:SigmaField},cellx::AbstractArray{<:VectorValue{2}})
  cache,  = c
  setsize!(cache,size(cellx))
  y = cache.array
  r = f.object.r

  map!(x -> TensorValue{2,3}(-r*sin( rem2pi(x[1],RoundDown) )*cos(x[2]),
                             -r*cos( rem2pi(x[1],RoundDown) )*sin(x[2]),
                              r*cos( rem2pi(x[1],RoundDown) )*cos(x[2]),
                             -r*sin( rem2pi(x[1],RoundDown) )*sin(x[2]),
                              0.0,
                              r*cos(x[2])),
                    y, cellx)

  return y

end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:SigmaField},x::VectorValue{2,T}) where {T}
  zero(TensorValue{2,3,T}) ## Jacobian is 3x2, recall we need to return transpose
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:SigmaField},x::VectorValue{2})
  y = cache
  r = f.object.r
  y = TensorValue{2,3}(-r*sin( rem2pi(x[1],RoundDown) )*cos(x[2]),
                       -r*cos( rem2pi(x[1],RoundDown) )*sin(x[2]),
                        r*cos( rem2pi(x[1],RoundDown) )*cos(x[2]),
                        -r*sin( rem2pi(x[1],RoundDown) )*sin(x[2]),
                        0.0,
                        r*cos(x[2]))
  return y
end


"""
SigmaMap
"""
struct SigmaMap{A} <: Map
  r::A # sphere radius
end

Sigma() = SigmaMap(r)


""" map 3D point on sphere -> latlon (2D) """
function Gridap.Arrays.return_cache(k::SigmaMap,cellx::AbstractArray{<:VectorValue{3,T}}) where {T}
  y = similar(cellx,VectorValue{2,T})
  return y #CachedArray(y)
end

function Gridap.Arrays.evaluate!(cache,f::SigmaMap,cellx::AbstractArray{<:VectorValue{3}})
  # 3D point on sphere -> lat lon
  y = cache
  map!(x -> VectorValue( rem2pi(atan(x[2], x[1]),RoundDown),
                                 atan(x[3], sqrt(x[1]*x[1] + x[2]*x[2]) )),
                        y, cellx)
  return y
end

function Gridap.Arrays.return_cache(k::SigmaMap,x::VectorValue{3,T}) where {T}
  y = zero(VectorValue{2,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SigmaMap,x::VectorValue{3})
  # 3D point on sphere -> lat lon
  # println("this func")
  y = cache
  y = VectorValue(  rem2pi(atan(x[2], x[1]),RoundDown),
                           atan(x[3], sqrt(x[1]*x[1] + x[2]*x[2]) ))
  return y
end

""" map latlon -> 3D point on sphere """
function Gridap.Arrays.return_cache(k::SigmaMap,latlon::AbstractArray{<:VectorValue{2,T}}) where {T}
  y = similar(latlon,VectorValue{3,T})
  return y#CachedArray(y)
end

function Gridap.Arrays.evaluate!(cache,f::SigmaMap,latlon::AbstractArray{<:VectorValue{2}})
  # lat-lon -> 3D point on sphere
  # latlon = θ, ϕ
  # setsize!(cache,size(cellx))
  y = cache
  r = f.r
  map!(x -> VectorValue(r*cos( rem2pi(x[1],RoundDown) )*cos(x[2]),
                        r*sin( rem2pi(x[1],RoundDown) )*cos(x[2]),
                        r*sin(x[2])),
                        y, latlon)
  return y
end

function Gridap.Arrays.return_cache(k::SigmaMap,x::VectorValue{2,T}) where {T}
  y = zero(VectorValue{3,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SigmaMap,x::VectorValue{2})
  # lat-lon -> 3D point on sphere
  # latlon = θ, ϕ

  y = cache
  r = f.r
  y = VectorValue(r*cos( rem2pi(x[1],RoundDown) )*cos(x[2]),
                  r*sin( rem2pi(x[1],RoundDown) )*cos(x[2]),
                  r*sin(x[2]))
  return y
end
