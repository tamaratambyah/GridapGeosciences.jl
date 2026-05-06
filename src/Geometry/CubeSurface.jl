"""
NPANELS

The number of panels of the cubed sphere manifold is always six.
This constant is used throughout the geometry construction
"""

const NPANELS = 6


"""
CUBE_HALF_EDGE

The value of the half edge of the cube is π/4. That is, every panel is the square
[-π/4,π/4]^2.
Currently the forward maps are only supported for this formulation of the panels.
"""

const CUBE_HALF_EDGE = π/4

function generate_ptr(Dc,n)
  nvertices = 2^Dc
  ptr  = Vector{Int}(undef,n+1)
  ptr[1]=1
  for i=1:n
    ptr[i+1]=ptr[i]+nvertices
  end
  ptr
end

function _CCAM_panel_wise_node_ids()
  ## CCAM panel ordering
  Dc=2
  data = [ 1,2,3,4, 3,4,5,6,  2,7,4,6, 8,5,7,6, 1,8,2,7, 1,3,8,5 ]
  ptr = generate_ptr(Dc,NPANELS)
  Table(data,ptr)
end

function _CCAM_cube_nodes_3d()
  CUBE_HALF_EDGE.* [
    Point(1.0, -1.0, -1.0)  # node 1
    Point(1.0, 1.0, -1.0)   # node 2
    Point(1.0, -1.0, 1.0)   # node 3
    Point(1.0, 1.0, 1.0)    # node 4
    Point(-1.0, -1.0, 1.0)  # node 5
    Point(-1.0, 1.0, 1.0)   # node 6
    Point(-1.0, 1.0, -1.0)  # node 7
    Point(-1.0, -1.0, -1.0) # node 8
   ]
end


function coarse_cube_surface_3D()

  nodes_3d = _CCAM_cube_nodes_3d()
  cell_node_ids = _CCAM_panel_wise_node_ids()

  polytopes = fill(QUAD,NPANELS)
  cell_type = fill(1,NPANELS)
  reffes = LagrangianRefFE(Float64,QUAD,1)
  cell_reffes=[reffes]

  topo = UnstructuredGridTopology(nodes_3d,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
  labels = FaceLabeling(topo)

  cube_grid = Gridap.Geometry.UnstructuredGrid(nodes_3d,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

  panel_ids = collect(1:NPANELS)
  return cube_grid,topo,labels,panel_ids
end



function coarse_cube_model()
  cube_grid,topo,labels, = coarse_cube_surface_3D()
  model = UnstructuredDiscreteModel(cube_grid,topo,labels)
  return model
end



function get_nodes_from_coords(topo::UnstructuredGridTopology{Dc,Dp},
  coords_array::AbstractArray{<:Vector{<:VectorValue{D,T}}}) where {Dc,Dp,D,T}

  cell_node_ids = get_faces(topo,Dc,0)
  nodes = similar(coords_array, VectorValue{D,T}, num_vertices(topo))

  get_nodes_from_coords!(nodes,cell_node_ids,coords_array)

  return nodes
end

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
