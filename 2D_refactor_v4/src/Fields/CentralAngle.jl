
"""
CentralAngleMap

Maps 2D local Cartesian points on the reference panel (panel 1) to the central
angles
cellx = (̃x,̃y)
central angles = (̃α,̃̃β)

̃α = atan(̃x)
̃β = atan(̃y)
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

"""
InverseCentralAngleMap

Maps central angles (̃α,̃̃β) to 2D local Cartesian points on the reference panel (̃x,̃y)
̃x = tan ̃α
̃y = tan ̃β
"""
struct InverseCentralAngleMap <: Map
end

function Gridap.Arrays.return_cache(k::InverseCentralAngleMap,cangles)
  y = similar(cangles)
  return y# CachedArray(y)
end

function Gridap.Arrays.evaluate!(cache,f::InverseCentralAngleMap,cangles)
  # setsize!(cache,size(cellx))
  y = cache
  map!(x -> VectorValue(tan(x[1]),tan(x[2])),
                        y, cangles)
  return y
end
