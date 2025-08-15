struct Cartesian2SphereicalMap  <: Field
end

"""
The forward map goes 3D -> 2D.
  θ = atan(Y,X)
  ϕ = atan(Z,sqrt(X^2 + Y^2))
"""

function Gridap.Arrays.return_cache(f::Cartesian2SphereicalMap,cellx::AbstractArray{<:VectorValue{3}})
  y = similar(cellx,VectorValue{2,Float64})
  return y
end


function Gridap.Arrays.evaluate!(cache,f::Cartesian2SphereicalMap,cellx::AbstractArray{<:VectorValue{3}} )

  y = cache
  map!(x-> cartesian_to_latlon(x), y, cellx  )

  return y
end


function Gridap.Arrays.return_cache(f::Cartesian2SphereicalMap,x::VectorValue{3})
  y = zero(VectorValue{2,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Cartesian2SphereicalMap,x::VectorValue{3})
  y = cache
  y = cartesian_to_latlon(x)
  return y
end
