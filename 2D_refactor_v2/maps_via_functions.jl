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
include("panel_rotations.jl")
include(srcdir("CoordinateTransformations.jl"))
include("bump_panel1.jl")
include("cube_surface_1_cell_per_panel_2D.jl")


"""
standard spherical map from lat-lon to points on the sphere
"""
struct sigmaMap <: Map
  xyz2θϕ::Function
  θϕ2xyz::Function
end

function Gridap.Arrays.evaluate!(cache,f::sigmaMap,θϕ::VectorValue{2})
  # input = lat-lon
  # returns a 3D point on the sphere
  f.θϕ2xyz(θϕ)
end

function Gridap.Arrays.evaluate!(cache,f::sigmaMap,X_s::VectorValue{3})
  # input = 3D point on sphere
  # returns lat-lon

  f.xyz2θϕ(X_s)
end

θϕ = Point(π,0.0)
sigma = sigmaMap(xyz2θϕ,θϕ2xyz)
Xs = evaluate(sigma,θϕ)
evaluate(sigma,Xs)
@test evaluate(sigma,Xs) == θϕ



"""
map between the reference panel (panel 1) and panels of the cube (1-6)
  requires panel_id as input
by default, map panel 1 -> panel p
to map panel p -> panel 1, set inverse = true
"""
struct panelMap <: Map
  rotate_panel_p_to_1::Vector{Matrix{Float64}}
  rotate_panel_1_to_p::Vector{Matrix{Float64}}
end


function Gridap.Arrays.evaluate!(cache,f::panelMap,X::VectorValue{3},panel_id::Int,inverse::Bool=false)
  A = f.rotate_panel_1_to_p[panel_id] # rotate panel 1 -> panel p

  if inverse # rotate panel p -> panel 1
    A = f.rotate_panel_p_to_1[panel_id]
  end

  TensorValue(A)⋅X
end

X = Point(1,1,1)
panel_map = panelMap(rotate_panel_p_to_1,rotate_panel_1_to_p)
evaluate!(nothing,panel_map,X,1)
evaluate!(nothing,panel_map,X,1,true)




"""
equi-angular gnomonic map to lat-lon
  input is local cartesian coords on reference panel
  returns lat-lon
"""
struct gammaMap <: Map
end

# map local cartesian coords on panel 1 to lat-lon
function Gridap.Arrays.evaluate!(cache,f::gammaMap,X::VectorValue{2})
  x,y = X # local 2D Cartesian coordinates on the reference panel

  @assert x∈[-1.0,1.0] && y∈[-1.0,1.0]

  θ = atan(x)

  bt = (1.0 + x*x + y*y)^(0.5)
  ϕ = asin( y/bt   )

  Point(θ,ϕ)
end

x = Point(-1,1)
gamma = gammaMap()
evaluate(gamma,x)


"""
map 2D->3D Cartesian on reference panel
"""
struct panel1BumpMap <: Map
  A::TensorValue{2,3}
  B::TensorValue{3,2}
  b::VectorValue{3}
end

function Gridap.Arrays.evaluate!(cache,f::panel1BumpMap,x::VectorValue{2})
  f.B ⋅ x + f.b
end

function Gridap.Arrays.evaluate!(cache,f::panel1BumpMap,X::VectorValue{3})
  f.A⋅X
end

A,B,b = bump_matrics()
panel1_bump = panel1BumpMap(A,B,b)

evaluate(panel1_bump,Point(1,1))
evaluate(panel1_bump,Point(1,1,1))

x = Point(0.0,1.0)
evaluate!(nothing,panel1_bump,x)

evaluate!(nothing,panel_map,evaluate!(nothing,panel1_bump,x),4 )



#### make a master map that combined all required maps

cube_grid,topo,face_labels = cube_surface_1_cell_per_panel_2D()
nodes = get_node_coordinates(cube_grid)
cell_ids = get_cell_node_ids(cube_grid)
ncells = length(cell_ids)
cell_map = get_cell_map(cube_grid)

evaluate(cell_map[1],Point(1,1))

panel_map = panelMap(rotate_panel_p_to_1,rotate_panel_1_to_p)
panel1_bump = panel1BumpMap( bump_matrics()...)

geo_map = (panel_map ∘panel1_bump)∘cell_map

master_map = Fill(Gridap.CellData.GenericField(geo_map),ncells)


evaluate(master_map[1],Point(1,1))
