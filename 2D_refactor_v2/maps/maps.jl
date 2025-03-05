## CachedArray, setsize!(cache,x) allocates? Ask Alberto

################################################################################
#####
"""
PanelRotationMap

Applies the roation provided by mats to x

`mats` can be either a single TensorValue, or an array of TensorValues
if `panel_id` is not provided, assume `mats` is the correct panel matrix

This map can also be used as the inverse, given mats^{-1} is the provided input
"""
struct PanelRotationMap{A} <: Map
  mats::A
end

function Gridap.Arrays.return_cache(f::PanelRotationMap,cellx::AbstractArray{<:VectorValue},panel_id::Int64)
  A = f.mats[panel_id]
  x = first(cellx)
  T = typeof(A⋅x)
  y = similar(cellx,T)
  return y # CachedArray(y)
end

function Gridap.Arrays.evaluate!(cache,f::PanelRotationMap,cellx::AbstractArray{<:VectorValue},panel_id::Int64)
  # setsize!(cache,size(cellx))
  y = cache
  A = f.mats[panel_id]
  map!(x -> A⋅x, y, cellx)
  return y
end

"""
if no panel id is provided, assume mats is the correct panel matrix
"""
function Gridap.Arrays.return_cache(f::PanelRotationMap,cellx::AbstractArray{<:VectorValue})
  A = f.mats
  x = first(cellx)
  T = typeof(A⋅x)
  y = similar(cellx,T)
  return y # CachedArray(y)
end

function Gridap.Arrays.evaluate!(cache,f::PanelRotationMap,cellx::AbstractArray{<:VectorValue})
  # setsize!(cache,size(cellx))
  y = cache
  A = f.mats
  map!(x -> A⋅x, y, cellx)
  return y
end

"""
if no panel id is provided, assume mats is the correct panel matrix
"""
function Gridap.Arrays.return_cache(f::PanelRotationMap,x::VectorValue)
  y = zero(x)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::PanelRotationMap,x::VectorValue)
  y = cache
  A = f.mats
  y = A .⋅ x
  return y
end

################################################################################
##### Bump map: panel 1 points 2D <-> 3D

"""
Panel1BumpMap

Map points on the reference panel (panel 1) from 2D <-> 3D.

The `A_bump`, `B_bump`, `b_bump` Tensor(Vector)Values are based on the rotation
of panel 1 in the C-CAM ordering

The dimension of cellx is used to determine what map to apply:
  if D == 3, bump 3D -> 2D
  if D == 2, bump 2D -> 3D

This map is the inverse of itself
"""

struct Panel1BumpMap{A,B,b} <: Map
  A_bump::A
  B_bump::B
  b_bump::b
end

function Gridap.Arrays.return_cache(f::Panel1BumpMap,
  cellx::AbstractArray{<:VectorValue{D}}) where {D}
  A = f.A_bump
  B = f.B_bump
  x = first(cellx)

  if D == 3 # D==3, -> y == 2 components; bump 3D -> 2D
    T = typeof(A⋅x)
  elseif D == 2 # D==2 -> y == 3 components;  bump 2D -> 3D
    T = typeof(B⋅x)
  end

  y = similar(cellx,T)
  return y# CachedArray(y)
end

function Gridap.Arrays.evaluate!(cache,f::Panel1BumpMap,
  cellx::AbstractArray{<:VectorValue{D}}) where {D}
  # setsize!(cache,size(cellx))
  y = cache
  A = f.A_bump
  B = f.B_bump
  b = f.b_bump

  if D == 3 # bump 3D -> 2D
    # y .= A.⋅cellx
    map!(x -> A⋅x, y, cellx)
  elseif D == 2 # bump 2D -> 3D
    # y .= B.⋅cellx .+ b
    map!(x -> B⋅x+b, y, cellx)
  end

  return y
end


function Gridap.Arrays.return_cache(f::Panel1BumpMap,x::VectorValue{D}) where {D}
  A = f.A_bump
  B = f.B_bump

  if D == 3 # D==3, -> y == 2 components; bump 3D -> 2D
    T = typeof(A⋅x)
  elseif D == 2 # D==2 -> y == 3 components;  bump 2D -> 3D
    T = typeof(B⋅x)
  end

  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Panel1BumpMap,x::VectorValue{D}) where {D}
  y = cache
  A = f.A_bump
  B = f.B_bump
  b = f.b_bump

  if D == 3 # bump 3D -> 2D
    y = A.⋅x
  elseif D == 2 # bump 2D -> 3D
    y = B.⋅x .+ b
  end

  return y
end


################################################################################
"""
GnomonicMap

Applies the Gnomonic maping from 2D local Cartesian points on the reference panel
(panel 1) to lat-lon on the reference panel.
See eq (32) in Giraldo et al. 2003, JCP, doi:10.1016/S0021-9991(03)00300-0

This implementation of the GnomonicMap is the composition of 2 maps,
  1. map 2D local Cartesian -> central angles on reference panel
  2. map central angles -> lat-lon on reference panel
"""
struct GnomonicMap <: Map
end

function Gridap.Arrays.return_cache(k::GnomonicMap,cellx)
  y = similar(cellx)
  return y# CachedArray(y)
end

function Gridap.Arrays.evaluate!(cache,f::GnomonicMap,cellx)
  # setsize!(cache,size(cellx))
  y = cache
  map!(x -> VectorValue(atan(x[1]),
                        asin( x[2]/( (1+x[1]^2 + x[2]^2)^(0.5) ) ) ),
                        y, cellx)
  return y
end


"""
CentralAngleMap

Maps 2D local Cartesian points on the reference panel (panel 1) to the central
angles
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
################################################################################
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
