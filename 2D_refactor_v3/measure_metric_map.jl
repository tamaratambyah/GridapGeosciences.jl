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


struct MeasureMult <: Map end

function Gridap.Arrays.return_cache(k::MeasureMult,bx::AbstractVector,gx::AbstractVector{<:TensorValue})
  g = first(gx)
  b = first(bx)
  T = typeof( sqrt(meas(g))*b  )
  _T = typeof(sqrt(meas(g)))
  y = similar(bx,T)
  z = similar(gx,_T)
  return y,z
end

function Gridap.Arrays.evaluate!(cache,k::MeasureMult,bx::AbstractVector,gx::AbstractVector{<:TensorValue})
  y, z = cache
  map!(x -> sqrt(meas(x)), z, gx)

  y .= bx .* z

  return y
end

_gx = fill( [TensorValue{2,2}(1.0,2.0,3.0,4.0) for i in 1:10], 10)
_bx = fill( [1.0 for i in 1:10], 10)

out = lazy_map(MeasureMult(), _bx, _gx)


cache = array_cache(out)
bm1() = lazy_collect(cache,out)
@benchmark bm1()
