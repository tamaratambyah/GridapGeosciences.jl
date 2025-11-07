struct Mapp{A}  <: Field
  n::A
end

function Gridap.Arrays.return_cache(f::Mapp,cellX::AbstractArray{<:VectorValue{3}})
  y = similar(cellX,VectorValue{3,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Mapp,cellX::AbstractArray{<:VectorValue{3}} )

  y = cache
  n = f.n

  y = map(cellX) do X
    γ = X⋅n - π/4# a
    return _intrusion(γ,X) + γ*n
  end

  return y
end

function Gridap.Arrays.return_cache(f::Mapp,X::VectorValue{3})
  T = typeof(x)
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Mapp,X::VectorValue{3})

  y = cache
  n = f.n

  γ = X⋅n - π/4# a
  y =  _intrusion(γ,X) + γ*n
  return y
end

_intrusion(γ::Float64,x) = x - γ*sqrt(3)*normal_vec(x)
