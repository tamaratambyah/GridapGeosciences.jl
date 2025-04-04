"""
Sigma
- SigmaMap
- SigmaField

The standard spherical mapping from lat-lon <-> 3D Cartesian on sphere surface.

The `latlon` variable is latlon = (θ, ϕ).
θ is the longitude, computed to be [-π,π]
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

struct SigmaField{A} <: Gridap.Fields.Field
  r::A # sphere radius
end


""" map 3D point on sphere -> latlon (2D) """
function Gridap.Arrays.return_cache(k::SigmaField,cellx::AbstractArray{<:VectorValue{3,T}}) where {T}
  y = similar(cellx,VectorValue{2,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SigmaField,cellx::AbstractArray{<:VectorValue{3}})
  # 3D point on sphere -> lat lon
  y = cache
  map!(x -> VectorValue( rem2pi(atan(x[2], x[1]),RoundNearest),  sin(x[3])), y, cellx)
  return y
end

function Gridap.Arrays.return_cache(k::SigmaField,x::VectorValue{3,T}) where {T}
  y = zero(VectorValue{2,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SigmaField,x::VectorValue{3})
  # 3D point on sphere -> lat lon
  y = cache
  y = VectorValue( rem2pi(atan(x[2], x[1]),RoundNearest),
                   sin(x[3]))
  return y
end

""" map latlon -> 3D point on sphere """
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
  map!(x -> VectorValue(r*cos(x[1])*cos(x[2]),
                        r*sin(x[1])*cos(x[2]),
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
  y = VectorValue(r*cos(x[1])*cos(x[2]),
                  r*sin(x[1])*cos(x[2]),
                  r*sin(x[2]))
  return y
end


## gradient: is just ....
# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
# Note  TensorValue{2,3}(0,4,1,0,5,1) == [0 1 5
#                                         4 0 1]
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:SigmaField},
  cellx::AbstractArray{<:VectorValue{2,T}}) where {T}
  _T = typeof(TensorValue{2,3,T})
  y = similar(cellx,T,size(cellx))
  CachedArray(y)
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:SigmaField},cellx::AbstractArray{<:VectorValue{2}})
  cache,  = c
  setsize!(cache,size(cellx))
  y = cache.array
  r = f.object.r

  map!(x -> TensorValue{2,3}(-r*sin(x[1])*cos(x[2]), -r*sin(x[2])*cos(x[1]),
                              r*cos(x[1])*cos(x[2]), -r*sin(x[2])*sin(x[1]),
                              0,                      r*cos(x[2])),
                    y, cellx)

  return y

end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:SigmaField},x::VectorValue{2,T}) where {T}
  zero(TensorValue{2,3,T}) ## Jacobian is 3x2, recall we need to return transpose
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:SigmaField},x::VectorValue{2})
  y = cache
  r = f.object.r
  y = TensorValue{2,3}(-r*sin(x[1])*cos(x[2]), -r*sin(x[2])*cos(x[1]),
                        r*cos(x[1])*cos(x[2]), -r*sin(x[2])*sin(x[1]),
                        0,                      r*cos(x[2]))
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
  map!(x -> VectorValue( rem2pi(atan(x[2], x[1]),RoundNearest),  sin(x[3])), y, cellx)
  return y
end

function Gridap.Arrays.return_cache(k::SigmaMap,x::VectorValue{3,T}) where {T}
  y = zero(VectorValue{2,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SigmaMap,x::VectorValue{3})
  # 3D point on sphere -> lat lon
  y = cache
  y = VectorValue( rem2pi(atan(x[2], x[1]),RoundNearest),
                   sin(x[3]))
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
  map!(x -> VectorValue(r*cos(x[1])*cos(x[2]),
                        r*sin(x[1])*cos(x[2]),
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
  y = VectorValue(r*cos(x[1])*cos(x[2]),
                  r*sin(x[1])*cos(x[2]),
                  r*sin(x[2]))
  return y
end
