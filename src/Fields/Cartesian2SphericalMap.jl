struct Cartesian2SphereicalMap <: Field
end

"""
The forward map goes 3D -> 2D.
  θ = atan(Y,X)
  ϕ = atan(Z,sqrt(X^2 + Y^2))

  Note, this map must be applied cellwise to obtain correct latlon.
  This is because if any(X > 0) and any(Y < 0) within the cell, take the θ angle
  from the negative side.
  Applying the map this way means that we account for the position of the cell within
  the whole sphere to obtain the correct longitude.
  Thus, the single point version of the map is not available
    ** update, have to make it available for vtk (unsure why)
"""

function Gridap.Arrays.return_cache(f::Cartesian2SphereicalMap,cellx::AbstractArray{<:VectorValue{3}})
  out = similar(cellx,VectorValue{2,Float64})
  x = similar(cellx,Float64)
  y = similar(cellx,Float64)
  z = similar(cellx,Float64)
  r = similar(cellx,Float64)
  θ = similar(cellx,Float64)
  ϕ = similar(cellx,Float64)
  return out,x,y,z,r,θ,ϕ
end


function Gridap.Arrays.evaluate!(cache,f::Cartesian2SphereicalMap,cellx::AbstractArray{<:VectorValue{3}} )
  # println("cell map")
  out,x,y,z,r,θ,ϕ = cache
  # map!(x-> xyz2θϕ(x,p), y, cellx  )

  # map!(q->q[1], x,cellx)
  # map!(q->q[2], y,cellx)
  # map!(q->q[3], z,cellx)
  # map!(x-> sqrt(x[1]^2 + x[2]^2 + x[3]^2) ,  r, cellx)

  # map!((x,y)->rem2pi(atan(y, x),RoundDown) ,  θ, x,y)
  # map!((z,r)-> asin(z/r), ϕ, z,r)

  x = map(x->x[1],cellx)
  y = map(x->x[2],cellx)
  z = map(x->x[3],cellx)
  r = sqrt.(x.^2 + y.^2 + z.^2)
  θ = rem2pi.(atan.(y, x),RoundDown)
  ϕ = asin.(z./r)

  # if there are negative ys and positive xs
  if any(y .< 0 ) && any(x .> 0)
    # map!((x,y)-> 2*π + rem2pi(atan(y, x),RoundUp) ,  θ, x,y)
    θ = 2*π .+  rem2pi.(atan.(y, x),RoundUp)
  end

  # map!((x,y)->VectorValue(x,y)  ,out, θ,ϕ)
  out = map(θ,ϕ) do θ,ϕ
    VectorValue(θ,ϕ)
  end

  return out
end


function Gridap.Arrays.return_cache(f::Cartesian2SphereicalMap,x::VectorValue{3})
  y = zero(VectorValue{2,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Cartesian2SphereicalMap,x::VectorValue{3})
  # println("single point")
  y = cache
  y = xyz2θϕ(x)
  return y
end
