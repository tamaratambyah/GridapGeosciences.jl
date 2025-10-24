struct Map2ExtrudedCube{A,B,C}  <: Field
  p::A
  Apanel::B
  b::C # z vector (0,0,1)
end


function Gridap.Arrays.return_cache(f::Map2ExtrudedPanel,cellx::AbstractArray{<:VectorValue{3}})
  x = first(cellx)
  T = typeof(x)
  y = similar(cellx,T)
  γ = similar(cellx,Float64)
  return y,γ
end


function Gridap.Arrays.evaluate!(cache,f::Map2ExtrudedPanel,cellx::AbstractArray{<:VectorValue{3}} )
  p = f.p
  y,γ = cache
  map!(x -> map_points_to_surface(p,x), y, cellx)
  return y
end


function Gridap.Arrays.return_cache(f::Map2ExtrudedPanel,x::VectorValue{3})
  T = typeof(x)
  y = zero(T)
  γ = 0.0
  return y,γ
end

function Gridap.Arrays.evaluate!(cache,f::Map2ExtrudedPanel,x::VectorValue{3})
  p = f.p
  Apanel = f.Apanel
  b = f.b
  y,γ = cache
  γ = extrusion_variable(p,x)
  y = _intrusion(γ,x) # surface points
  return Apanel ⋅ y + γ*b
end
