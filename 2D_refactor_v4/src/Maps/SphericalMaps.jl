"""
GnomonicMap

Applies the Gnomonic maping from central angles to lat-lon on the reference panel.
See eq (32) in Giraldo et al. 2003, JCP, doi:10.1016/S0021-9991(03)00300-0
"""
struct GnomonicMap <: Map
end

function Gridap.Arrays.return_cache(k::GnomonicMap,cangles)
  y = similar(cangles)
  return y# CachedArray(y)
end

function Gridap.Arrays.evaluate!(cache,f::GnomonicMap,cangles)
  # setsize!(cache,size(cellx))
  y = cache
  map!(x -> VectorValue(x[1],
                        asin( ( tan(x[2]) ) / ( (1+ (tan(x[1]))^2 + (tan(x[2]))^2 )^(0.5) ) )
                        ),
                        y, cangles)
  return y
end


"""
SigmaMap

The standard spherical mapping from lat-lon <-> 3D Cartesian on sphere surface.

The `latlon` variable is latlon = (θ, ϕ).
θ is the longitude, computed to be [-π,π]
ϕ is the latitude, computed to be [-π/2, π/2]

`r` is the radius of the sphere, which is assumed to be constant

This map is the inverse of itself, dispatch based on number of components
  i.e. cellx = VectorValue{3}, latlon = VectorValue{2}
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
  # setsize!(cache,size(cellx))
  y = cache
  map!(x -> VectorValue( rem2pi(atan(x[2], x[1]),RoundNearest),  sin(x[3])), y, cellx)
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

################################################################################
"""
CentralAngleMap

Maps 2D local Cartesian points on the reference panel (panel 1) to the central
angles
cellx = (̃x,̃y)
central angles = (̃α,̃̃β)

̃α = atan(̃x)
̃β = atan(̃y)
"""
struct CentralAngleMap <: Map
end

function Gridap.Arrays.return_cache(k::CentralAngleMap,cellx)
  y = similar(cellx)
  return y# CachedArray(y)
end

function Gridap.Arrays.evaluate!(cache,f::CentralAngleMap,cellx)
  # setsize!(cache,size(cellx))
  y = cache
  map!(x -> VectorValue(atan(x[1]),atan(x[2])),
                        y, cellx)
  return y
end

"""
InverseCentralAngleMap

Maps central angles (̃α,̃̃β) to 2D local Cartesian points on the reference panel (̃x,̃y)
̃x = tan ̃α
̃y = tan ̃β
"""
struct InverseCentralAngleMap <: Map
end

function Gridap.Arrays.return_cache(k::InverseCentralAngleMap,cangles)
  y = similar(cangles)
  return y# CachedArray(y)
end

function Gridap.Arrays.evaluate!(cache,f::InverseCentralAngleMap,cangles)
  # setsize!(cache,size(cellx))
  y = cache
  map!(x -> VectorValue(tan(x[1]),tan(x[2])),
                        y, cangles)
  return y
end
