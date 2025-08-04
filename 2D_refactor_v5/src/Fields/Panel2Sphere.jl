
## panels 1,2,3
struct Panel123Sphere{A} <: Gridap.Fields.Field
  r::A # sphere radius
end


""" map x,y -> 3D point on sphere for xy local in panel123"""
function Gridap.Arrays.return_cache(k::Panel123Sphere,cellx::AbstractArray{<:VectorValue{2,T}}) where {T}
  y = similar(cellx,VectorValue{3,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Panel123Sphere,cellx::AbstractArray{<:VectorValue{2}})
  y = cache
  r = f.r

  map!(x -> VectorValue(r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(0.5) ) * 1.0,
                        r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(0.5) ) * tan(x[1]),
                        r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(0.5) ) * tan(x[2])),
                        y, cellx)
  return y
end

function Gridap.Arrays.return_cache(k::Panel123Sphere,x::VectorValue{2,T}) where {T}
  y = zero(VectorValue{3,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Panel123Sphere,x::VectorValue{2})
  # lat-lon -> 3D point on sphere
  # latlon = θ, ϕ

  y = cache
  r = f.r
  y = VectorValue(r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(0.5) ) * 1.0,
                  r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(0.5) ) * tan(x[1]),
                  r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(0.5) ) * tan(x[2]))
  return y
end


"""
gradient of local x,y in panel1,2,3 -> 3D point on sphere
  J = r/ρ^3 ( -tan(x)sec^2(y)       -tan(y)sec^(y)
              sec^2(x)sec^(y)       -tan(x)tan(y)sec^(y)
              -tan(y)tan(x)sec^(x)  sec^2(y)sec^2(x)
       )
where ρ^3 = (1 + tan^2(x) + tan^2(y) )^(3/2)

# Gridap convention dictates we return the transpose (https://github.com/gridap/Gridap.jl/issues/822)
# Note  TensorValue{2,3}(0,4,1,0,5,1) == [0 1 5
                                          4 0 1]

The transpose is:
  JT = r/ρ^3 ( -tan(x)sec^2(y)  sec^2(x)sec^(y)       -tan(y)tan(x)sec^(x)
               -tan(y)sec^(y)  -tan(x)tan(y)sec^(y)   sec^2(y)sec^2(x) )

To write as a TensorValue:
  TensorValue{2,3}  r/ρ^3 ( -tan(x)sec^2(y),  -tan(y)sec^(y), sec^2(x)sec^(y), -tan(x)tan(y)sec^(y), -tan(y)tan(x)sec^(x), sec^2(y)sec^2(x) )
"""
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:Panel123Sphere},
  cellx::AbstractArray{<:VectorValue{2,T}}) where {T}
  _T = typeof(TensorValue{2,3,T})
  y = similar(cellx,_T,size(cellx))
  CachedArray(y)
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:Panel123Sphere},cellx::AbstractArray{<:VectorValue{2}})
  cache,  = c
  setsize!(cache,size(cellx))
  y = cache.array
  r = f.object.r


  map!(x -> TensorValue{2,3}( r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*(sec(x[1]))^2,
                              r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*(sec(x[2]))^2,
                              r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *  (sec(x[1]))^2*(sec(x[2]))^2,
                              r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*tan(x[2])*(sec(x[2]))^2,
                              r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*tan(x[1])*(sec(x[1]))^2,
                              r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *   (sec(x[1]))^2*(sec(x[2]))^2
  ),
                    y, cellx)

  return y

end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:Panel123Sphere},x::VectorValue{2,T}) where {T}
  zero(TensorValue{2,3,T}) ## Jacobian is 3x2, recall we need to return transpose
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:Panel123Sphere},x::VectorValue{2})
  y = cache

  r = f.object.r
  y = TensorValue{2,3}( r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*(sec(x[1]))^2,
                        r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*(sec(x[2]))^2,
                        r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *  (sec(x[1]))^2*(sec(x[2]))^2,
                        r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*tan(x[2])*(sec(x[2]))^2,
                        r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*tan(x[1])*(sec(x[1]))^2,
                        r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *   (sec(x[1]))^2*(sec(x[2]))^2  )
  return y
end





## panels 4, 5, 6
struct Panel456Sphere{A} <: Gridap.Fields.Field
  r::A # sphere radius
end


""" map x,y -> 3D point on sphere for xy local in panel123"""
function Gridap.Arrays.return_cache(k::Panel456Sphere,cellx::AbstractArray{<:VectorValue{2,T}}) where {T}
  y = similar(cellx,VectorValue{3,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Panel456Sphere,cellx::AbstractArray{<:VectorValue{2}})
  y = cache
  r = f.r

  map!(x -> VectorValue(r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(0.5) ) * 1.0,
                        r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(0.5) ) * tan(x[2]),
                        r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(0.5) ) * tan(x[1])),
                        y, cellx)
  return y
end

function Gridap.Arrays.return_cache(k::Panel456Sphere,x::VectorValue{2,T}) where {T}
  y = zero(VectorValue{3,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Panel456Sphere,x::VectorValue{2})
  # lat-lon -> 3D point on sphere
  # latlon = θ, ϕ

  y = cache
  r = f.r
  y = VectorValue(r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(0.5) ) * 1.0,
                  r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(0.5) ) * tan(x[2]),
                  r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(0.5) ) * tan(x[1]))
  return y
end



""" jacobian"""
function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:Panel456Sphere},
  cellx::AbstractArray{<:VectorValue{2,T}}) where {T}
  _T = typeof(TensorValue{2,3,T})
  y = similar(cellx,_T,size(cellx))
  CachedArray(y)
end

function Gridap.Arrays.evaluate!(c,f::FieldGradient{1,<:Panel456Sphere},cellx::AbstractArray{<:VectorValue{2}})
  cache,  = c
  setsize!(cache,size(cellx))
  y = cache.array
  r = f.object.r

  map!(x -> TensorValue{2,3}( r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*(sec(x[1]))^2,
                              r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*(sec(x[2]))^2,
                              r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*tan(x[2])*(sec(x[1]))^2,
                              r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *  (sec(x[1]))^2*(sec(x[2]))^2,
                              r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *  (sec(x[1]))^2*(sec(x[2]))^2,
                              r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*tan(x[1])*(sec(x[2]))^2
  ),
                    y, cellx)

  return y

end


function Gridap.Arrays.return_cache(cache,f::FieldGradient{1,<:Panel456Sphere},x::VectorValue{2,T}) where {T}
  zero(TensorValue{2,3,T}) ## Jacobian is 3x2, recall we need to return transpose
end

function Gridap.Arrays.evaluate!(cache,f::FieldGradient{1,<:Panel456Sphere},x::VectorValue{2})
  y = cache
  r = f.object.r
  y =  TensorValue{2,3}(  r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*(sec(x[1]))^2,
                          r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*(sec(x[2]))^2,
                          r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[1])*tan(x[2])*(sec(x[1]))^2,
                          r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *  (sec(x[1]))^2*(sec(x[2]))^2,
                          r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) *  (sec(x[1]))^2*(sec(x[2]))^2,
                          r/( (1 + (tan(x[1]))^2 + (tan(x[2]))^2 )^(1.5) ) * -1.0*tan(x[2])*tan(x[1])*(sec(x[2]))^2  )
  return y
end
