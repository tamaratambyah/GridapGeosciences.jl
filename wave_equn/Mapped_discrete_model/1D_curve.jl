"""
Take a 1D point and map it to a 2D point
The 1D point can be thought of as the parametrising coordinate.
i.e. a circle can be parameterised by a single coordinate, θ, as [sin(θ),cos(θ)]
"""

using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Test

include("feproblem.jl")

map(x) = VectorValue(x...,4*x...)
# map(x) = VectorValue(cos(x[1]...),sin(x[1]...))
evaluate(map,Point(1.0))

p = 1
u_ex(x) = 2*x[1]

domain = (0,4,)
partition = (2)
model = CartesianDiscreteModel(domain,partition)
grid = get_grid(model)

uh = fe_problem(model,p,u_ex)
writevtk(Triangulation(model),datadir("Mapped")*"/1D_curve_og",cellfields=["u"=>uh],append=false)

# new_grid = MappedGrid(get_grid(model),map)



cmaps = get_cell_map(grid)
evaluate(cmaps[1],Point(1))

cell_node_ids = get_cell_node_ids(grid)
old_nodes = get_node_coordinates(grid)



dims = num_dims(grid) + 1
T = VectorValue{dims,Float64} # eltype(old_nodes)

node_coordinates = Vector{T}(undef,length(old_nodes))
c_coor = get_cell_coordinates(grid)./1

N = length(c_coor)
map_c_coor  = fill(Vector{T}(undef,length(c_coor)),N)
for i = 1:N
  # map_c_coor[i] = Base.map( (x)->evaluate(map,x), c_coor[i])
  map_c_coor[i] = Base.map( (x)->map(x), c_coor[i])
end
println(map_c_coor)
println(c_coor)
g = lazy_map(  (x)-> map(x), c_coor)


# Gridap.Arrays.lazy_map(map,c_coor)
Gridap.Geometry._cell_vector_to_dof_vector!(node_coordinates,cell_node_ids,map_c_coor)
node_coords = collect(node_coordinates)

model_map=get_cell_map(grid)
# geo_map=lazy_map(∘,map,model_map)

import Gridap.Geometry.Grid


struct MyMappedGrid{Dc,Dp,M,L} <: Grid{Dc,Dp}
  grid::Grid{Dc,Dp}
  phys_map::M # New map in the physical space
  node_coords::L
  function MyMappedGrid(grid::Grid{Dc,Dp},phys_map) where {Dc,Dp}

    cell_node_ids = get_cell_node_ids(grid)
    old_nodes = get_node_coordinates(grid)

    dims = num_dims(grid) + 1
    T = VectorValue{dims,Float64} # eltype(old_nodes)

    node_coordinates = Vector{T}(undef,length(old_nodes))
    c_coor = get_cell_coordinates(grid)./1

    N = length(c_coor)
    map_c_coor  = fill(Vector{T}(undef,length(c_coor)),N)
    for i = 1:N
      map_c_coor[i] = Base.map( (x)->evaluate(map,x), c_coor[i])
    end

    Gridap.Geometry._cell_vector_to_dof_vector!(node_coordinates,cell_node_ids,map_c_coor)
    node_coords = collect(node_coordinates)

    new{Dc,Dp,typeof(phys_map),typeof(node_coords)}(grid,phys_map,node_coords)
  end
end


Gridap.Geometry.get_node_coordinates(grid::MyMappedGrid) = grid.node_coords
Gridap.Geometry.get_cell_node_ids(grid::MyMappedGrid) = get_cell_node_ids(grid.grid)
Gridap.Geometry.get_reffes(grid::MyMappedGrid) = get_reffes(grid.grid)
Gridap.Geometry.get_cell_type(grid::MyMappedGrid) = get_cell_type(grid.grid)




struct MyMappedDiscreteModel{Dc,Dp} <: DiscreteModel{Dc,Dp}
  model::DiscreteModel{Dc,Dp}
  mapped_grid
  function MyMappedDiscreteModel(model::DiscreteModel{Dc,Dp},node_coords) where {Dc,Dp}
    new{Dc,Dp}(model,node_coords)
  end
end

Gridap.Geometry.get_grid(model::MyMappedDiscreteModel) = model.mapped_grid
Gridap.Geometry.get_cell_map(model::MyMappedDiscreteModel) = get_cell_map(model.mapped_grid)
Gridap.Geometry.get_grid_topology(model::MyMappedDiscreteModel) = Gridap.Geometry.get_grid_topology(model.model)
Gridap.Geometry.get_face_labeling(model::MyMappedDiscreteModel) = get_face_labeling(model.model)

mapp_grid = MyMappedGrid(grid,map)

new_model = MyMappedDiscreteModel(model,mapp_grid)


function Gridap.Visualization._vtkpoints(trian)
  D = num_point_dims(trian)
  if D == 1
    D += 1
  end
  println(D)
  x = get_node_coordinates(trian)
  xflat = collect(x)
  reshape(reinterpret(Float64,xflat),(D,length(x)))
end

writevtk(new_model,datadir("Mapped")*"/1D_curve",append=false)


uh_map = fe_problem(new_model,p,u_ex)
writevtk(Triangulation(new_model),datadir("Mapped")*"/1D_curve_mapped",cellfields=["u_map"=>uh_map],append=false)
