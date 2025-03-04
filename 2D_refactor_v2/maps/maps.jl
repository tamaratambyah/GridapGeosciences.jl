################################################################################
##### Rotation map
struct PanelRotationMap{A} <: Map # rotation panel p using mats[p]
  mats::A
end

function Gridap.Arrays.return_cache(f::PanelRotationMap,cellx,panel_id::Int64)
  A = f.mats[panel_id]
  x = first(cellx)
  T = typeof(A⋅x)
  y = similar(cellx,T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::PanelRotationMap,cellx,panel_id::Int64)
  y = cache
  A = f.mats[panel_id]
  map!(x -> A⋅x, y, cellx)
  return y
end



################################################################################
##### Bump map: panel 1 points 2D <-> 3D
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
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Panel1BumpMap,
  cellx::AbstractArray{<:VectorValue{D}}) where {D}
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



################################################################################
##### Gnomonic  mapping 2D local Cartesian on panel 1 -> latlon
struct GnomonicMap <: Map
end

function Gridap.Arrays.return_cache(k::GnomonicMap,cellx)
  y = similar(cellx)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::GnomonicMap,cellx)
  y = cache
  map!(x -> VectorValue(atan(x[1]),
                        asin( x[2]/( (1+x[1]^2 + x[2]^2)^(0.5) ) ) ),
                        y, cellx)
  return y
end


################################################################################
##### sigma standard spherical mapping from latlon <-> 3D Cartesian
struct SigmaMap{A} <: Map
  r::A # sphere radius
end

# map 3D point on sphere -> latlon (2D)
function Gridap.Arrays.return_cache(k::SigmaMap,cellx::AbstractArray{<:VectorValue{3,T}}) where {T}
  y = similar(cellx,VectorValue{2,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SigmaMap,cellx::AbstractArray{<:VectorValue{3}})
  # 3D point on sphere -> lat lon
  y = cache
  map!(x -> VectorValue( rem2pi(atan(x[2], x[1]),RoundNearest),  sin(x[3])), y, cellx)
  return y
end

# map latlon -> 3D point on sphere
function Gridap.Arrays.return_cache(k::SigmaMap,latlon::AbstractArray{<:VectorValue{2,T}}) where {T}
  y = similar(latlon,VectorValue{3,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SigmaMap,latlon::AbstractArray{<:VectorValue{2}})
  # lat-lon -> 3D point on sphere
  # latlon = θ, ϕ
  y = cache
  r = f.r
  map!(x -> VectorValue(r*cos(x[1])*cos(x[2]),
                        r*sin(x[1])*cos(x[2]),
                        r*sin(x[2])),
                        y, latlon)
  return y
end
