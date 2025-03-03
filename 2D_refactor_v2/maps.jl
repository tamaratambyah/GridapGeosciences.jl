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
using StaticArrays
include("panel_rotations.jl")
include("bump_panel1.jl")
# include("cube_surface_1_cell_per_panel_2D.jl")

function lazy_collect(cache,arr)
  for i in eachindex(arr)
    getindex!(cache, arr, i)
  end
end


"""
standard spherical map from lat-lon to points on the sphere
"""
struct SigmaMap{T} <: Map
  r::T
end

function Gridap.Arrays.return_cache(k::SigmaMap,θϕ::VectorValue{2})
  T = eltype(θϕ)
  D = length(θϕ)
  y = zero(VectorValue{D+1,T})
  return y
end


function Gridap.Arrays.evaluate!(cache,f::SigmaMap,θϕ::VectorValue{2})
  # input = lat-lon
  # returns a 3D point on the sphere
  θ,ϕ = θϕ
  x,y,z = cache
  r = f.r

  x = r.*cos.(θ).*cos.(ϕ)
  y = r.*sin.(θ).*cos.(ϕ)
  z = r.*sin.(ϕ)

  VectorValue(x,y,z)
end

struct InvSigmaMap <: Map
end


function Gridap.Arrays.return_cache(k::InvSigmaMap,Xs::VectorValue{3})
  T = eltype(Xs)
  D = length(Xs)
  y = zero(VectorValue{D-1,T})
  return y
end


function Gridap.Arrays.evaluate!(cache,f::InvSigmaMap,Xs::VectorValue{3})
  # input = lat-lon
  # returns a 3D point on the sphere
  θ,ϕ = cache
  x,y,z = Xs
  θ = atan.(y, x)
  ϕ = asin.(z)
  VectorValue(θ,ϕ)

end


# map lat-lon -> sphere
a = 1.0
r = a/sqrt(3)
θϕs = [Point(π,0.0) for i in 1:1000]
sigma = lazy_map(SigmaMap(r), θϕs)
cache = array_cache(sigma)
bm1() = lazy_collect(cache,sigma)
@benchmark bm1()

# map sphere -> lat-lon
Xs = [Point(-1/sqrt(3),0.0,0.0) for i in 1:1000]
inv_sigma = lazy_map(InvSigmaMap(), Xs)
cache = array_cache(inv_sigma)
bm2() = lazy_collect(cache,inv_sigma)
@benchmark bm2()

# test composition, should be inverse
Master = Gridap.Arrays.Operation(SigmaMap(r))(InvSigmaMap())
z = lazy_map(Master,Xs)
cache = array_cache(z)
bm3() = lazy_collect(cache,z)
@benchmark bm3()

Master = Gridap.Arrays.Operation(InvSigmaMap())(SigmaMap(r))
z = lazy_map(Master,θϕs)
cache = array_cache(z)
bm4() = lazy_collect(cache,z)
@benchmark bm4()

###############################################################################


"""
map between the reference panel (panel 1) and panels of the cube (1-6)
  requires panel_id as input
by default, map panel 1 -> panel p
to map panel p -> panel 1, set inverse = true
"""
const rotate_panel_p_to_1, rotate_panel_1_to_p = panel_rotations()

struct PanelMap <: Map end # rotate panel 1 -> panel p

function Gridap.Arrays.return_cache(f::PanelMap,X::VectorValue{3},panel_id::Int)
  y = VectorValue(get_array(X)) # data from X into y
  A =  TensorValue( rotate_panel_1_to_p[panel_id] ) # rotate panel 1 -> panel p
  return y, A
end

function Gridap.Arrays.evaluate!(cache,f::PanelMap,X::VectorValue{3},panel_id::Int)
  y, A = cache
  y = A⋅ X
end


struct InvPanelMap <: Map end # rotate panel p -> panel 1

function Gridap.Arrays.return_cache(f::InvPanelMap,X::VectorValue{3},panel_id::Int)
  y = VectorValue(get_array(X)) # copy data from X into y
  A = TensorValue( rotate_panel_p_to_1[panel_id] ) # rotate panel p -> panel 1
 return y, A
end

function Gridap.Arrays.evaluate!(cache,f::InvPanelMap,X::VectorValue{3},panel_id::Int)
  y, A = cache
  y = A⋅ X
end

panel_id = [rand(1:6) for i in 1:1000]
Xs =  [Point(1.0,1.0,1.0) for i in 1:1000]
panelMap = lazy_map(PanelMap(), Xs, panel_id)
cache = array_cache(panelMap)

bm1() = lazy_collect(cache,panelMap)
@benchmark bm1()


invPanelMap = lazy_map(InvPanelMap(), Xs, panel_id)
cache = array_cache(invPanelMap)

bm2() = lazy_collect(cache,invPanelMap)
@benchmark bm2()


# compose - should return X
Master = Gridap.Arrays.Operation(InvPanelMap())(PanelMap())
z = lazy_map(Master,Xs,panel_id)
cache = array_cache(z)
bm3() = lazy_collect(cache,z)
@benchmark bm3()



"""
equi-angular gnomonic map to lat-lon
  input is local cartesian coords on reference panel
  returns lat-lon
"""
struct gammaMap <: Map
end

# map local cartesian coords on panel 1 to lat-lon
function Gridap.Arrays.return_cache(f::gammaMap,X::VectorValue{2})
  y = zero(X)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::gammaMap,X::VectorValue{2})
  θϕ = cache
  x,y = X # local 2D Cartesian coordinates on the reference panel

  @assert (-1.0<=x<=1.0)  && ( -1.0<=y<=1.0 )

  θ = atan.(x) # theta

  bt = (1.0 + x.*x + y.*y).^(0.5)
  ϕ = asin( y./bt   ) # phi

  θϕ = VectorValue(θ,ϕ)

  return θϕ
end

xs = [Point(0.0,1.0) for i in 1:1000]
gamma = lazy_map(gammaMap(),xs)
cache = array_cache(gamma)
bm1() = lazy_collect(cache,gamma)
@benchmark bm1()





"""
map 2D->3D Cartesian on reference panel
"""
const A_bump,B_bump,b_bump = bump_matrics()

struct panel1BumpMap <: Map
end

function Gridap.Arrays.return_cache(f::panel1BumpMap,X::VectorValue{2})
  y = zero(X)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::panel1BumpMap,x::VectorValue{2})
  y = cache
  y = B_bump⋅x + b_bump
  return y
end

function Gridap.Arrays.return_cache(f::panel1BumpMap,X::VectorValue{3})
  y = zero(X)
  return y
end

function Gridap.Arrays.evaluate!(cache,f::panel1BumpMap,X::VectorValue{3})
  y = cache
  y = A_bump⋅X
  return y
end

xs = [Point(0.0,1.0) for i in 1:1000]
bump = lazy_map(panel1BumpMap(),xs)
cache = array_cache(bump)
bm1() = lazy_collect(cache,bump)
@benchmark bm1()

Xs = [Point(1.0,1.0,1.0) for i in 1:1000]
bump = lazy_map(panel1BumpMap(),Xs)
cache = array_cache(bump)
bm2() = lazy_collect(cache,bump)
@benchmark bm2()


# compose - should return x=1
Master = Gridap.Arrays.Operation(panel1BumpMap())(panel1BumpMap())
z = lazy_map(Master,Xs)
cache = array_cache(z)
bm3() = lazy_collect(cache,z)
@benchmark bm3()

z = lazy_map(Master,xs)
cache = array_cache(z)
bm4() = lazy_collect(cache,z)
@benchmark bm4()
