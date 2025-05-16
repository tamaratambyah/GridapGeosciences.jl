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




struct SigmaField1{A} <: Gridap.Fields.Field
  r::A # sphere radius
end



""" map latlon -> 3D point on sphere, σ: (θ,ϕ) → (X_s, Y_s, Z_s) """
function Gridap.Arrays.return_cache(k::SigmaField1,latlon::AbstractArray{<:VectorValue{2,T}}) where {T}
  y = similar(latlon,VectorValue{3,T})
  return y#CachedArray(y)
end

function Gridap.Arrays.evaluate!(cache,f::SigmaField1,latlon::AbstractArray{<:VectorValue{2}})
  # lat-lon -> 3D point on sphere
  # latlon = θ, ϕ
  # setsize!(cache,size(cellx))
  y = cache
  r = f.r
  map!(x -> VectorValue(r*cos( (x[1]) )*cos(x[2]),
                        r*sin( (x[1]) )*cos(x[2]),
                        r*sin(x[2])),
                        y, latlon)
  return y
end

function Gridap.Arrays.return_cache(k::SigmaField1,x::VectorValue{2,T}) where {T}
  y = zero(VectorValue{3,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SigmaField1,x::VectorValue{2})
  # lat-lon -> 3D point on sphere
  # latlon = θ, ϕ

  y = cache
  r = f.r
  y = VectorValue(r*cos( (x[1]) )*cos(x[2]),
                  r*sin( (x[1]) )*cos(x[2]),
                  r*sin(x[2]))
  return y
end





################################################################################
#### Sigma 2
################################################################################

struct SigmaField2{A} <: Gridap.Fields.Field
  r::A # sphere radius
end



""" map latlon -> 3D point on sphere, σ: (θ,ϕ) → (X_s, Y_s, Z_s) """
function Gridap.Arrays.return_cache(k::SigmaField2,latlon::AbstractArray{<:VectorValue{2,T}}) where {T}
  y = similar(latlon,VectorValue{3,T})
  return y#CachedArray(y)
end

function Gridap.Arrays.evaluate!(cache,f::SigmaField2,latlon::AbstractArray{<:VectorValue{2}})
  y = cache
  r = f.r
  map!(x -> VectorValue(r*cos( (x[1]) )*sin(x[2]),
                        r*sin( (x[1]) )*sin(x[2]),
                        r*cos(x[2])),
                        y, latlon)
  return y
end

function Gridap.Arrays.return_cache(k::SigmaField2,x::VectorValue{2,T}) where {T}
  y = zero(VectorValue{3,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SigmaField2,x::VectorValue{2})
  # lat-lon -> 3D point on sphere
  # latlon = θ, ϕ

  y = cache
  r = f.r
    y = VectorValue(r*cos( (x[1]) )*sin(x[2]),
                    r*sin( (x[1]) )*sin(x[2]),
                    r*cos(x[2]))
  return y
end
