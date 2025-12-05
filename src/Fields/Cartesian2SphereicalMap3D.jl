struct Cartesian2SphereicalMap3D <: Field
end

"""
The forward map goes 3D -> 3D.
  θ = atan(Y,X)
  ϕ = atan(Z,sqrt(X^2 + Y^2))
  r = sqrt(X^2 + Y^2 + Z^2)

  For future consideration: remove Cartesian2SphereicalMap, and obtain θϕ as
  [θ   = [1 0 0   [θ
   ϕ]     0 1 0]   ϕ
                   r]
"""

function Gridap.Arrays.return_cache(f::Cartesian2SphereicalMap3D,cellx::AbstractArray{<:VectorValue{3}})
  out = similar(cellx,VectorValue{3,Float64})
  return out
end


function Gridap.Arrays.evaluate!(cache,f::Cartesian2SphereicalMap3D,cellx::AbstractArray{<:VectorValue{3}} )
  # println("cell map")
  out = cache

  x = map(x->x[1],cellx)
  y = map(x->x[2],cellx)
  z = map(x->x[3],cellx)

  ## This hack is because sometimes at higher refinements, the y is slightly negative/positive
  ## when it really should be zero. It is numerical error. So to overcome, set
  ## the abs(y) <1e-16 (machine eps) to be zero
  idx = abs.(y) .<1e-16
  if any(idx)
    y[idx].= 0.0
  end

  r = sqrt.(x.^2 + y.^2 + z.^2)
  θ = rem2pi.(atan.(y, x),RoundDown)
  ϕ = asin.(z./r)

  # if there are negative ys and positive xs
  if any(y .< 0 ) && any(x .> 0)
    θ = 2*π .+  rem2pi.(atan.(y, x),RoundUp)
  end



  # map!((x,y)->VectorValue(x,y)  ,out, θ,ϕ)
  out = map(θ,ϕ,r) do θ,ϕ,r
    VectorValue(θ,ϕ,r)
  end

  return out
end


function Gridap.Arrays.return_cache(f::Cartesian2SphereicalMap3D,x::VectorValue{3})
  y = zero(VectorValue{3,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Cartesian2SphereicalMap3D,x::VectorValue{3})
  # println("single point")
  y = cache
  y = xyz2θϕr(x)
  return y
end
