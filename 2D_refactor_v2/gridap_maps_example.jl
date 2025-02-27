using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools

function lazy_collect(cache,arr)
  for i in eachindex(arr)
    getindex!(cache, arr, i)
  end
end

xs = [[1.0,4.0] for i in 1:2000]

f(x) = x.*x
g(x) = 2.0.*x


y = f.(xs)
@allocated f.(xs)

bm0() = f.(xs)
@benchmark bm0()

y = map(f, xs)
@allocated map(f, xs)
@allocated map(x-> f(x),xs)
bm1() = map(f,xs)
@benchmark bm1()

@allocated lazy_map(f, xs)
@allocated collect(lazy_map(f, xs))

y = lazy_map(f,xs)
cache = array_cache(y)
bm2() = lazy_collect(cache,y)
@benchmark bm2()
@allocated bm2()

struct SquareMap <: Map
end

function Gridap.Arrays.return_cache(k::SquareMap,x::Vector)
  y = similar(x)
  return y
end

function Gridap.Arrays.evaluate!(cache,k::SquareMap,x::Vector)
  y = cache
  y .= x.*x
  return y
end

y = lazy_map(SquareMap(),xs)
cache = array_cache(y)

bm3() = lazy_collect(cache,y)
@benchmark bm3()


struct Times2Map <: Map
end

function Gridap.Arrays.return_cache(k::Times2Map,x::Vector)
  y = similar(x)
  return y
end

function Gridap.Arrays.evaluate!(cache,k::Times2Map,x::Vector)
  y = cache
  y .= 2.0.*x
  return y
end

z = lazy_map(Times2Map(),xs)
cache = array_cache(y)
bm4() = lazy_collect(cache,z)

@benchmark bm4()


# master
h = g∘f
master0() = h.(xs)
@benchmark master0()

y = lazy_map(h,xs)
cache = array_cache(y)

master1() = lazy_collect(cache,y)
@benchmark master1()

Master = Gridap.Arrays.Operation(g)(f)
z = lazy_map(Master,xs)
cache = array_cache(z)

master3() = lazy_collect(cache,z)

@benchmark master3()

# composed functions

Master = Gridap.Arrays.Operation(Times2Map())(SquareMap())
z = lazy_map(Master,xs)
cache = array_cache(z)
master4() = lazy_collect(cache,z)

@benchmark master4()
