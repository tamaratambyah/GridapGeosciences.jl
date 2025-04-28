"""
Gnomonic
- GnomonicMap
- GnomonicField
- InvGnomonicMap
- InvGnomonicField

Applies the Gnomonic maping from central angles to lat-lon on the reference panel.
- See eq (32) in Giraldo et al. 2003, JCP, doi:10.1016/S0021-9991(03)00300-0
- Note, `cangle` variable is central_angle = (α, β), where α,β ∈ [-π/4, π/4]^2

The same functionality is provided for a GridapMap and GridapField
- the GridapField is used to map cell coordinates and cell maps
    * the gradient is implemented
- the GridapMap is generally not used, but is provided regardless for testing
purposes
"""

"""
GnomonicField

γ: (α,β) → (θ,ϕ)
map central angles (2D) → latlon (2D)
"""

struct GnomonicField <: Gridap.Fields.Field
end

function Gridap.Arrays.return_cache(k::GnomonicField,cellcangles::AbstractArray{<:VectorValue{2}})
  y = similar(cellcangles)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::GnomonicField,cellcangles::AbstractArray{<:VectorValue{2}})
  # setsize!(cache,size(cellx))
  y = cache
  map!(x -> VectorValue(x[1],
                        asin( ( tan(x[2]) ) / ( (1+ (tan(x[1]))^2 + (tan(x[2]))^2 )^(0.5) ) )
                        ),
                        y, cellcangles)
  return y
end



function Gridap.Arrays.return_cache(k::GnomonicField,cangle::VectorValue{2})
  y = zero(VectorValue{2,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::GnomonicField,cangle::VectorValue{2})
  y = cache
  y =  VectorValue(cangle[1],
                  asin( ( tan(cangle[2]) ) / ( (1+ (tan(cangle[1]))^2 + (tan(cangle[2]))^2 )^(0.5) ) )
                        )
  return y
end

"""
gradient of γ: (α,β) → (θ,ϕ) is:
  J = [ 1 0
        F G ]
where
F = tan(α)*sec(α)*tan(β)/ρ^2
G = (sec(β))^2 sec(α)/ρ^2
ρ^2 = 1 + (tanα)^2 + (tanβ)^2

# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
# Note  TensorValue{2,3}(0,4,1,0,5,1) == [0 1 5
                                          4 0 1]
"""

function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:GnomonicField},cangles::AbstractArray{<:VectorValue{2}})
  _T = typeof(TensorValue{2,2,Float64})
  y = similar(cangles,_T,size(cangles))
  CachedArray(y)
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:GnomonicField},cangles::AbstractArray{<:VectorValue{2}})
  cache,  = c
  setsize!(cache,size(cangles))
  y = cache.array

  map!(x -> TensorValue{2,2}( 1, (-1/( 1 + (tan(x[1]))^2 + (tan(x[2]))^2) )*tan(x[1])*sec(x[1])*tan(x[2]),
                              0, ( 1/( 1 + (tan(x[1]))^2 + (tan(x[2]))^2) )*(sec(x[2]))^2*sec(x[1])
                            ),
                    y, cangles)

  return y
end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:GnomonicField},x::VectorValue{2})
  zero(TensorValue{2,2,Float64})
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:GnomonicField},x::VectorValue{2})
  y = cache
  y = TensorValue{2,2}( 1, (-1/( 1 + (tan(x[1]))^2 + (tan(x[2]))^2) )*tan(x[1])*sec(x[1])*tan(x[2]),
                        0, ( 1/( 1 + (tan(x[1]))^2 + (tan(x[2]))^2) )*(sec(x[2]))^2*sec(x[1])
                                     )
  return y
end

# function Gridap.Arrays.return_cache(cache,f::FieldGradient{2,<:GnomonicField},x::VectorValue{2})
#   zero(TensorValue{2,2,Float64})
# end

# function Gridap.Arrays.evaluate!(cache,f::FieldGradient{2,<:GnomonicField},x::VectorValue{2})
#   y = cache
#   y = zero(TensorValue{2,2,Float64})
#   return y
# end


"""
InvGnomonicField

γ^{-1}: (θ,ϕ) → (α,β)
map latlon (2D) → central angles (2D)
"""

struct InvGnomonicField <: Gridap.Fields.Field
end

function Gridap.Arrays.return_cache(k::InvGnomonicField,latlons::AbstractArray{<:VectorValue{2}})
  y = similar(latlons)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::InvGnomonicField,latlons::AbstractArray{<:VectorValue{2}})
  # setsize!(cache,size(cellx))
  y = cache
  map!(x -> VectorValue(latlon[1],
                        atan( tan(latlon[2]), cos(latlon[1]) )
                        ),
                        y, latlons)
  return y
end



function Gridap.Arrays.return_cache(k::InvGnomonicField,latlon::VectorValue{2})
  y = zero(VectorValue{2,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::InvGnomonicField,latlon::VectorValue{2})
  y = cache
  y =  VectorValue(latlon[1],
                   atan( tan(latlon[2]), cos(latlon[1]) )
                        )
  return y
end


"""
gradient of γ^{-1}: (θ,ϕ) → (α,β) is:
  J = [ 1 0
        F G ]
where
F = (cosθ*tanϕ*tanθ)/( (tanϕ)^2 + (cosθ)^2 )
G = cosθ/( ( (tanϕ)^2 + (cosθ)^2  )*(cosϕ)^2   )

# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
# Note  TensorValue{2,3}(0,4,1,0,5,1) == [0 1 5
                                          4 0 1]
"""


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:InvGnomonicField},latlons::AbstractArray{<:VectorValue{2}})
  _T = typeof(TensorValue{2,2,Float64})
  y = similar(latlons,_T,size(latlons))
  CachedArray(y)
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:InvGnomonicField},latlons::AbstractArray{<:VectorValue{2}})
  cache,  = c
  setsize!(cache,size(latlons))
  y = cache.array

  map!(x -> TensorValue{2,2}( 1, ( cos(x[1])*tan(x[2])*tan(x[1]) )/( (tan(x[2]))^2 + (cos(x[1]))^2 ),
                              0, ( cos(x[1]) )/( (cos(x[2]))^2*( (tan(x[2]))^2 + (cos(x[1]))^2 ) )
                              ),
                    y, latlons)

  return y
end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:InvGnomonicField},x::VectorValue{2})
  zero(TensorValue{2,2,Float64})
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:InvGnomonicField},x::VectorValue{2})
  y = cache
  y = TensorValue{2,2}( 1, ( cos(x[1])*tan(x[2])*tan(x[1]) )/( (tan(x[2]))^2 + (cos(x[1]))^2 ),
                        0, ( cos(x[1]) )/( (cos(x[2]))^2*( (tan(x[2]))^2 + (cos(x[1]))^2 ) )
                                     )
  return y
end





"""
GnomonicMap
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
