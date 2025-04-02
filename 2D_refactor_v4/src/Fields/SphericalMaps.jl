"""
GnomonicMap

Applies the Gnomonic maping from central angles to lat-lon on the reference panel.
See eq (32) in Giraldo et al. 2003, JCP, doi:10.1016/S0021-9991(03)00300-0

Note, `cangle` variable is central_angle = (α, β), where α,β ∈ [-π/4, π/4]^2
"""
struct GnomonicMap <: Map
end

function Gridap.Arrays.return_cache(k::GnomonicMap,cellcangles::AbstractArray{<:VectorValue{2}})
  y = similar(cellcangles)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::GnomonicMap,cellcangles::AbstractArray{<:VectorValue{2}})
  # setsize!(cache,size(cellx))
  y = cache
  map!(x -> VectorValue(x[1],
                        asin( ( tan(x[2]) ) / ( (1+ (tan(x[1]))^2 + (tan(x[2]))^2 )^(0.5) ) )
                        ),
                        y, cellcangles)
  return y
end



function Gridap.Arrays.return_cache(k::GnomonicMap,cangle::VectorValue{2})
  y = zero(VectorValue{2,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::GnomonicMap,cangle::VectorValue{2})
  y = cache
  y =  VectorValue(cangle[1],
                  asin( ( tan(cangle[2]) ) / ( (1+ (tan(cangle[1]))^2 + (tan(cangle[2]))^2 )^(0.5) ) )
                        )
  return y
end



"""
InvGnomonicMap

Applies the inverse Gnomonic maping from lat-lon to central angles on the
reference panel.
See eq (32) in Giraldo et al. 2003, JCP, doi:10.1016/S0021-9991(03)00300-0

Note, `latlon` variable is latlon = (θ, ϕ).
θ is the longitude, computed to be [-π,π]
ϕ is the latitude, computed to be [-π/2, π/2]
"""

struct InvGnomonicMap <: Map
end

function Gridap.Arrays.return_cache(k::InvGnomonicMap,latlons::AbstractArray{<:VectorValue{2}})
  y = similar(latlons)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::InvGnomonicMap,latlons::AbstractArray{<:VectorValue{2}})
  y = cache
  map!(x -> VectorValue(x[1],
                        atan( tan(x[2]), cos(x[1]))
                        ),
                        y, latlons)
  return y
end


function Gridap.Arrays.return_cache(k::InvGnomonicMap,latlon::VectorValue{2,T}) where {T}
  y = zero(VectorValue{2,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::InvGnomonicMap,latlon::VectorValue{2})
  y = cache
  y = VectorValue(latlon[1],
                  atan( tan(latlon[2]), cos(latlon[1]))
                        )
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
