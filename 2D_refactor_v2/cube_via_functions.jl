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

include("maps/panel_rotations.jl")
include("maps/bump_panel1.jl")
include("cube_topo/cube_surface_1_cell_per_panel.jl")
include("maps/panel_ids_from_refinement_v2.jl")


dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

function lazy_collect(cache,arr)
  for i in eachindex(arr)
    getindex!(cache, arr, i)
  end
end


function make_grid(topo::UnstructuredGridTopology{Dc,Dp},nodes::AbstractArray) where {Dc,Dp}
  # Dc = num_cell_dims(topo)
  cell_reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
  cell_node_ids = get_faces(topo,Dc,0)
  cell_type = get_cell_type(topo)

  Gridap.Geometry.UnstructuredGrid(nodes,cell_node_ids,cell_reffes,cell_type,OrientationStyle(topo))
end



"""
map between the reference panel (panel 1) and panels of the cube (1-6)
  requires panel_id as input
by default, map panel 1 -> panel p
to map panel p -> panel 1, set inverse = true
"""

function panel_map(panel_id::Int64,inverse::Bool=false)
  rotate_panel_p_to_1, rotate_panel_1_to_p = panel_rotations()

  A = rotate_panel_1_to_p[panel_id] # rotate panel 1 -> panel p

  if inverse # rotate panel p -> panel 1
    A = rotate_panel_p_to_1[panel_id]
  end

  function tmp(X)
    TensorValue(A)⋅X
  end
end


"""
map 2D<->3D Cartesian on reference panel
"""
function panel1_2D_to_3D(x::VectorValue{2,Float64})
  A,B,b = bump_matrics()
  B⋅x .+ b
end

function panel1_3D_to_2D(x::VectorValue{3,Float64})
  A,B,b = bump_matrics()
  A⋅x
end

# array of 6 maps
paneli3D_2_panel12D = [Gridap.Arrays.Operation(panel1_3D_to_2D)(panel_map(panel_id,true)) for panel_id in 1:6]
panel12D_2_paneli3D = [Gridap.Arrays.Operation(panel_map(panel_id,false))(panel1_2D_to_3D) for panel_id in 1:6]

cube_model_3D = UnstructuredDiscreteModel(cube_surface_1_cell_per_panel()...)
Dc = num_cell_dims(cube_model_3D)
Dp = num_cell_dims(cube_model_3D)

# nodes = get_node_coordinates(cube_model_3D)

# panel_id = 4
# z = lazy_map(paneli3D_2_panel12D[panel_id],nodes)
# cache = array_cache(z)
# master() = lazy_collect(cache,z)
# @benchmark master()

# mapped_nodes = collect(z)
# z = lazy_map(panel12D_2_paneli3D[panel_id],mapped_nodes)
# cache = array_cache(z)
# master() = lazy_collect(cache,z)
# @benchmark master()

# _nodes = collect(z)



###############################################################################
# test coarsest model
model = cube_model_3D
panel_ids = get_panel_ids(model)
cmaps = get_cell_map(model)
nodes = get_node_coordinates(model)

# map array of maps
f = [paneli3D_2_panel12D[i] for i in panel_ids ]
finv = [panel12D_2_paneli3D[i] for i in panel_ids ]
master_map = lazy_map(∘,f,cmaps)

# test for node 4 = (1,1,1) in cells 1,2,3
N4 = nodes[4]
cell, ref_point = 1, Point(1,1)
x = evaluate(master_map[cell],ref_point)
@test evaluate(finv[cell],x) == N4

cell, ref_point = 2, Point(1,0)
x = evaluate(master_map[cell],ref_point)
@test evaluate(finv[cell],x) == N4

cell, ref_point = 3, Point(0,0)
x = evaluate(master_map[cell],ref_point)
@test evaluate(finv[cell],x) == N4

### map all cells
mapped_nodes = Vector{VectorValue{Dc+1,Float64}}(undef,length(nodes))
ref_cell_coords = get_cell_ref_coordinates(get_grid(model))
ncells = num_cells(model)
cell_ids = get_cell_node_ids(model)


for i in 1:ncells

  phys_cell_coords_2D = map(x->evaluate(master_map[i],x),ref_cell_coords[i]  )

  # check the points are actually on the 2D ref panel
  for i in 1:length(phys_cell_coords_2D)
    x,y = phys_cell_coords_2D[i]
    # println("x = $x; y = $y")
    @assert (-1.0<=x<=1.0)  && ( -1.0<=y<=1.0 )
  end


  phys_cell_coords_3D = map(x->evaluate(finv[i],x),phys_cell_coords_2D  )


  @test phys_cell_coords_3D == nodes[cell_ids[i]]

  mapped_nodes[cell_ids[i]] .= 1.0.*phys_cell_coords_3D

end


mapped_grid = make_grid(get_grid_topology(cube_model_3D),mapped_nodes)
writevtk(mapped_grid,dir*"/cube_model",append=false)




################################################################################
# test a refined model
model = Gridap.Adaptivity.refine(cube_model_3D)
panel_ids = get_panel_ids(model)
cmaps = get_cell_map(model)
nodes = get_node_coordinates(model)

# map array of maps
f = [paneli3D_2_panel12D[i] for i in panel_ids ]
finv = [panel12D_2_paneli3D[i] for i in panel_ids ]
master_map = lazy_map(∘,f,cmaps)


### map all cells
mapped_nodes = Vector{VectorValue{Dc+1,Float64}}(undef,length(nodes))
ref_cell_coords = get_cell_ref_coordinates(get_grid(model))
ncells = num_cells(model)
cell_ids = get_cell_node_ids(model)

for i in 1:ncells

  phys_cell_coords_2D = map(x->evaluate(master_map[i],x),ref_cell_coords[i]  )

  # check the points are actually on the 2D ref panel
  for i in 1:length(phys_cell_coords_2D)
    x,y = phys_cell_coords_2D[i]
    # println("x = $x; y = $y")
    @assert (-1.0<=x<=1.0)  && ( -1.0<=y<=1.0 )
  end


  phys_cell_coords_3D = map(x->evaluate(finv[i],x),phys_cell_coords_2D  )


  @test phys_cell_coords_3D == nodes[cell_ids[i]]

  mapped_nodes[cell_ids[i]] .= 1.0.*phys_cell_coords_3D

end


mapped_grid = make_grid(get_grid_topology(model),mapped_nodes)
writevtk(mapped_grid,dir*"/ref_cube_model",append=false)
