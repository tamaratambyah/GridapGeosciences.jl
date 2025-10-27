struct Mapp  <: Field
end

function Gridap.Arrays.return_cache(f::Mapp,cellαβγ::AbstractArray{<:VectorValue{3}})
  y = similar(cellαβγ,VectorValue{3,Float64})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Mapp,cellαβγ::AbstractArray{<:VectorValue{3}} )

  y = cache
  map!(x-> remove_extrusion(x),
      y, cellαβγ  )
  return y
end

function Gridap.Arrays.return_cache(f::Mapp,αβγ::VectorValue{3})
  T = typeof(x)
  y = zero(T)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::Mapp,αβγ::VectorValue{3})

  y = cache
  y = remove_extrusion(x) # surface points
  return y
end

function remove_extrusion(αβγ::VectorValue{3}; a=π/4)
  α,β,γ = αβγ

  α1 = sign(α)*(abs(α) - 1.0)
  β1 = sign(β)*(abs(β) - 1.0)

  if γ == 0.0
    return αβγ
  else
    return Point(α1,β1,γ)
  end


end
