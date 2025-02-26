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

function sigma(θϕ<:VectorValue{2})
  θϕ2xyz(θϕ)
end

function sigma(X_s<:VectorValue{3})
  xyz2θϕ(X_s)
end

"""
map between the reference panel (panel 1) and panels of the cube (1-6)
  requires panel_id as input
by default, map panel 1 -> panel p
to map panel p -> panel 1, set inverse = true
"""

function panel_map(panel_id::Int,X<:VectorValue{3},inverse::Bool=false)
  A = rotate_panel_1_to_p[panel_id] # rotate panel 1 -> panel p

  if inverse # rotate panel p -> panel 1
    A = rotate_panel_p_to_1[panel_id]
  end

  TensorValue(A)⋅X
end


"""
equi-angular gnomonic map to lat-lon
  input is local cartesian coords on reference panel
  returns lat-lon
"""
function gamma(X<:VectorValue{2})
  x,y = X # local 2D Cartesian coordinates on the reference panel

  @assert x∈[-1.0,1.0] && y∈[-1.0,1.0]

  θ = atan(x)

  bt = (1.0 + x*x + y*y)^(0.5)
  ϕ = asin( y/bt   )

  Point(θ,ϕ)
end


"""
map 2D->3D Cartesian on reference panel
"""
function bump_panel1(x<:VectorValue{2})
  A,B,b = bump_matrics()
  B⋅x + b
end

function bump_panel1(x<:VectorValue{3})
  A,B,b = bump_matrics()
  A⋅x
end


#### make a master map that combined all required maps

cube_grid,topo,face_labels = cube_surface_1_cell_per_panel_2D()
nodes = get_node_coordinates(cube_grid)
cell_ids = get_cell_node_ids(cube_grid)
ncells = length(cell_ids)
cell_map = get_cell_map(cube_grid)

evaluate(cell_map[5],Point(0,0))

pids = collect(1:6)


_geo_map(p) = x ->  panel_map(p,bump_panel1(x))
geo_map = Fill(Gridap.CellData.GenericField(_geo_map(pids)),ncells)
master_map = lazy_map(∘,geo_map,cell_map)


evaluate(master_map[1],Point(1,1))
