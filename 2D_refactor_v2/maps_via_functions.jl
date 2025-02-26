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
include("cube_surface_1_cell_per_panel.jl")

"""
standard spherical map from lat-lon to points on the sphere
"""

function sigma(θϕ::VectorValue{2,Float64})
  θϕ2xyz(θϕ)
end

function sigma(X_s::VectorValue{3,Float64})
  xyz2θϕ(X_s)
end

"""
map between the reference panel (panel 1) and panels of the cube (1-6)
  requires panel_id as input
by default, map panel 1 -> panel p
to map panel p -> panel 1, set inverse = true
"""

function panel_map(panel_id::Int64,X::VectorValue{3,Float64},inverse::Bool=false)
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
function gamma(X::VectorValue{2,Float64})
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
function bump_panel1(x::VectorValue{2,Float64})
  A,B,b = bump_matrics()
  B⋅x + b
end

function bump_panel1(x::VectorValue{3,Float64})
  A,B,b = bump_matrics()
  A⋅x
end


#### make a master map that combined all required maps
pids = collect(1:6)
cube_model_2D = UnstructuredDiscreteModel(cube_surface_1_cell_per_panel_2D()...)
cube_model_3D = UnstructuredDiscreteModel(cube_surface_1_cell_per_panel()...)

ref_coord = Point(1,1) # point on reference FE (0,0)->(1,1)
cell_id = 4 # for testing

cmaps_2D = get_cell_map(cube_model_2D)
ncells = num_cells(cube_model_2D)

evaluate(cmaps_2D[4],ref_coord) # zero (as expected)

_geo_map(p) = x ->  panel_map(p,bump_panel1(x))
f = lazy_map(p->_geo_map(p),pids)

# make the all cell maps the same as cell 1
# this will return a point on panel 1
_cell_map = Fill(Gridap.CellData.GenericField(cmaps_2D[1]),ncells)
master_map = lazy_map(∘,f,_cell_map)

ref_panel_point = evaluate(_cell_map[cell_id],ref_coord) ### evaluate at ref FE point
cube_panel_point = evaluate(f[cell_id],ref_panel_point) ### evaluate at point on panel 1

evaluate(master_map[cell_id],ref_coord)

# test against the 3D cells which are 'correct'
cmaps_3D = get_cell_map(cube_model_3D)
ref_cell_coords = get_cell_ref_coordinates(get_grid(cube_model_2D))

for i in 1:ncells
  phys_cell_coords_2D = map(x->evaluate(master_map[i],x),ref_cell_coords[i]  )
  phys_cell_coords_3D = map(x->evaluate(cmaps_3D[i],x),ref_cell_coords[i]  )

  @test phys_cell_coords_2D == phys_cell_coords_3D
end
