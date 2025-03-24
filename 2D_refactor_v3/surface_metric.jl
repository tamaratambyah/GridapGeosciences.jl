"""
Consider surface metric of the form:
g = [E F
     F G]
"""

function factor(x)
  α,β = x
   r^2/( (1 + (tan(α))^2 + (tan(β))^2 )^2 * (cos(α))^2 * (cos(β))^2 )
end

function E(x)
  α,β = x
  factor(x)*( 1 + (tan(α))^2 )
end

function F(x)
  α,β = x
  -1.0*factor(x)*( tan(α)*tan(β)  )
end

function G(x)
  α,β = x
  factor(x)*( 1 + (tan(β))^2 )
end

sqrt_det_func(x) = sqrt(  E(x)*G(x) - (F(x))^2  )



using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.Fields
using Gridap.TensorValues
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools
include("helpers.jl")


r = 1*sqrt(3)
struct SurfaceMetric{A} <: Map
  r::A
end

function Gridap.Arrays.return_cache(k::SurfaceMetric,x::VectorValue{D,T}) where {D,T}
  y = zero(TensorValue{D,D,T})
  return y
end

function Gridap.Arrays.evaluate!(cache,f::SurfaceMetric,x::VectorValue{D,T}) where {D,T}

  y = cache
  r = f.r
  α,β = x

  factor = r^2/( (1 + (tan(α))^2 + (tan(β))^2 )^2 * (cos(α))^2 * (cos(β))^2 )

  E = factor*( 1 + (tan(α))^2 )
  F = -1.0*factor*( tan(α)*tan(β)  )
  G = factor*( 1 + (tan(β))^2 )

  y = TensorValue{D,D}(E,F,F,G)
  return y
end

x = [Point(1.0,1.0) for i in 1:6]
z = lazy_map(SurfaceMetric(r),x)
cache = array_cache(z)
bm() = lazy_collect(cache,z)
@benchmark bm()
