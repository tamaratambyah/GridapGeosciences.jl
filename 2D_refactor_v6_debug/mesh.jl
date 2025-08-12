using DrWatson
using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Helpers
using LinearAlgebra

dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


global RADIUS = 1.0

include("forward_map.jl")
include("inverse_map.jl")
include("MatMultField.jl")
include("../2D_refactor_v6/panel_ids_from_refinement.jl")

function generate_ptr(n)
  nvertices = 4
  ptr  = Vector{Int}(undef,n+1)
  ptr[1]=1
  for i=1:n
    ptr[i+1]=ptr[i]+nvertices
  end
  ptr
end


function coarse_cube_surface_3D(a::Real,npanels)

  nodes_3d = a.* [
    Point(1.0, -1.0, -1.0)  # node 1
    Point(1.0, 1.0, -1.0)   # node 2
    Point(1.0, -1.0, 1.0)   # node 3
    Point(1.0, 1.0, 1.0)    # node 4
    Point(-1.0, -1.0, 1.0)  # node 5
    Point(-1.0, 1.0, 1.0)   # node 6
    Point(-1.0, 1.0, -1.0)  # node 7
    Point(-1.0, -1.0, -1.0) # node 8
  ]

  ## CCAM panel ordering
  data = [ 1,2,3,4, 3,4,5,6,  2,7,4,6, 8,5,7,6, 1,8,2,7, 1,3,8,5 ]

  ptr = generate_ptr(npanels)
  cell_node_ids = Table(data,ptr)

  polytopes = fill(QUAD,npanels)
  cell_type = fill(1,npanels)
  reffes = LagrangianRefFE(Float64,QUAD,1)
  cell_reffes=[reffes]

  topo = UnstructuredGridTopology(nodes_3d,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
  labels = FaceLabeling(topo)

  cube_grid = Gridap.Geometry.UnstructuredGrid(nodes_3d,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

  panel_ids = collect(1:npanels)
  return cube_grid,topo,labels,panel_ids
end


npanels = 6
cube_grid,topo,labels, = coarse_cube_surface_3D(π/4,npanels)
cube_model = UnstructuredDiscreteModel(cube_grid,topo,labels)

cube_model = Gridap.Adaptivity.refine(cube_model)
cube_model = Gridap.Adaptivity.refine(cube_model)

writevtk(Triangulation(cube_model),dir*"/cube_mode",append=false)


n = Int(num_cells(cube_model)/npanels)
panel_ids = get_panel_ids(cube_model)

################################################################################
## make panel grid
cube_grid = get_grid(cube_model)
cube_topo = get_grid_topology(cube_model)

nodes = get_node_coordinates(cube_grid)
cmaps = get_cell_map(cube_grid)

swap = lazy_map(p->MatMultField( cube_to_αβ(p) ),panel_ids)

new_cmaps = lazy_map(∘,swap,cmaps)

ref_points = get_cell_ref_coordinates(cube_grid)
lazy_map(evaluate,new_cmaps,ref_points)./1

new_nodes = map(x->Point(x[2],x[3]),nodes)
new_grid = Geometry.UnstructuredGrid(new_nodes,get_cell_node_ids(cube_grid),get_reffes(cube_grid),get_cell_type(cube_grid),OrientationStyle(cube_grid),
                    nothing,new_cmaps)

new_topo = UnstructuredGridTopology(new_nodes,get_cell_node_ids(cube_grid),get_cell_type(cube_topo),get_polytopes(cube_topo),OrientationStyle(cube_topo))
new_labels = FaceLabeling(new_topo)
panel_model = UnstructuredDiscreteModel(new_grid,new_topo,new_labels)

################################################################################
## make ambient grid
function get_nodes_from_coords(grid::Grid{Dc,Dp},
  coords_array::AbstractArray{<:Vector{<:VectorValue{D,T}}}) where {Dc,Dp,D,T}

  cell_node_ids = get_cell_node_ids(grid)
  nodes = similar(coords_array, VectorValue{D,T}, num_nodes(grid))

  get_nodes_from_coords!(nodes,cell_node_ids,coords_array)

  return nodes
end

function get_nodes_from_coords!(nodes,cell_node_ids,
  coords_array::AbstractArray{<:Vector{<:VectorValue}})

  cache = array_cache(coords_array)

  for i in eachindex(coords_array)
    ids = cell_node_ids[i]
    nodes[ids] .= getindex!(cache, coords_array, i)
  end

end

panel_grid = get_grid(panel_model)
panel_topo = get_grid_topology(panel_model)

cmaps = get_cell_map(panel_grid)

panel_to_ambient_map = lazy_map(p->ForwardMap(p), panel_ids)

ambient_maps = lazy_map(∘,panel_to_ambient_map,cmaps)

ref_points = get_cell_ref_coordinates(panel_grid)
ambient_cell_coords = lazy_map(evaluate,ambient_maps,ref_points)./1

ambient_nodes = get_nodes_from_coords(panel_grid,ambient_cell_coords)

ambient_grid = Geometry.UnstructuredGrid(ambient_nodes,get_cell_node_ids(panel_grid),get_reffes(panel_grid),get_cell_type(panel_grid),OrientationStyle(panel_grid),
                    nothing,ambient_maps)

ambient_topo = UnstructuredGridTopology(ambient_nodes,get_cell_node_ids(panel_grid),get_cell_type(panel_topo),get_polytopes(panel_topo),OrientationStyle(panel_topo))
ambient_labels = FaceLabeling(panel_topo)

ambient_model = UnstructuredDiscreteModel(ambient_grid,ambient_topo,ambient_labels)

writevtk(Triangulation(ambient_model),dir*"/ambient_model",append=false)
