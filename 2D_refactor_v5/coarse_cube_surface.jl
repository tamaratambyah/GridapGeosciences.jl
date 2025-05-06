"""
coarse_cube_surface_2D -- builds the flat cube geometry, with correct
periodicity, with 1 cell per panel based on the CCAM ordering of panels. This
geometry contains 8 nodes. The following table information is required:

This cube has faces [-a,a] where the side length is 2a.
For the standard cube with faces [-1,1], set a = 1.
For the central angle mesh where cube faces are [-π/4,π/4]^2, set a = π/4.


panel/cell no | cell_node_ids
      1       |   [1 2 3 4]
      2       |   [3 4 5 6]
      3       |   [4 2 6 7]
      4       |   [6 7 5 8]
      5       |   [7 2 8 1]
      6       |   [8 1 5 3]

  node no     |   Point in 2D   |   Point in 3D
      1       |   a * (-1,-1)   |   a * (1,-1,-1)
      2       |   a * (1,-1)    |   a * (1,1,-1)
      3       |   a * (-1,1)    |   a * (1,-1,1)
      4       |   a * (1,1)     |   a * (1,1,1)
      5       |   a * (0,0)     |   a * (-1,-1,1)
      6       |   a * (0,0)     |   a * (-1,1,1)
      7       |   a * (0,0)     |   a * (-1,1,-1)
      8       |   a * (0,0)     |   a * (-1,-1,-1)

Note, the nodes are 2D. This is so that Dp = 2. Only the nodes on panel 1
(nodes 1-4) make sense in the nodes_2d array. The remaining nodes are set to
(0,0). The 3D nodes on panel p can be obtained by applying the bump 2D->3D map,
and then the rotation map to panel 1.
"""


function coarse_cube_surface_3D(a::Float64)
  npanels = 6

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
  data = [ 1,2,3,4, 3,4,5,6, 2,7,4,6, 8,5,7,6, 1,8,2,7, 1,3,8,5  ] # reorient + rotated

  ptr = generate_ptr(npanels)
  cell_node_ids = Table(data,ptr)

  polytopes = fill(QUAD,npanels)
  cell_type = fill(1,npanels)
  reffes = LagrangianRefFE(Float64,QUAD,1)
  cell_reffes=[reffes]

  topo = UnstructuredGridTopology(nodes_3d,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
  face_labels = FaceLabeling(topo)

  cube_grid = Gridap.Geometry.UnstructuredGrid(nodes_3d,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

  panel_ids = collect(1:6)
  return cube_grid,topo,face_labels,panel_ids
end


function get_cube_nodes_2D(cube_cell_coords_3D,panel_ids)
  coords_panel1_3D = lazy_map(Rp1PanelMap3D(), cube_cell_coords_3D, panel_ids)
  coords_panel1_2D = lazy_map(BumpMap(), coords_panel1_3D)
  return coords_panel1_2D
end

function cube_surface_2D(cube_grid_3D::Grid{Dc,Dp},
  topo_3D::GridTopology{Dc,Dp},panel_ids::Vector{Int}) where {Dc,Dp}

  cell_coords_3D = get_cell_coordinates(cube_grid_3D)
  cell_coords_panel1_2D = get_cube_nodes_2D(cell_coords_3D,panel_ids)

  polytopes = map(get_polytope, get_reffes(cube_grid_3D))

  nodes_2D = get_panel_1_nodes_from_coords(cube_grid_3D,cell_coords_panel1_2D,panel_ids)

  topo_2D = UnstructuredGridTopology(nodes_2D,
      get_cell_node_ids(cube_grid_3D),get_cell_type(cube_grid_3D),polytopes,Gridap.Geometry.NonOriented())
  face_labels_2D = FaceLabeling(topo_2D)
  cube_grid_2D = Gridap.Geometry.UnstructuredGrid(nodes_2D,
        get_cell_node_ids(cube_grid_3D),get_reffes(cube_grid_3D),get_cell_type(cube_grid_3D),Gridap.Geometry.NonOriented())

  @check (num_point_dims(cube_grid_2D) == num_point_dims(topo_2D) == Dp-1)
  @check (num_cell_dims(cube_grid_2D) == num_cell_dims(topo_2D) == Dc)

  return cube_grid_2D,topo_2D,face_labels_2D,panel_ids
end



function generate_ptr(n)
  nvertices = 4
  ptr  = Vector{Int}(undef,n+1)
  ptr[1]=1
  for i=1:n
    ptr[i+1]=ptr[i]+nvertices
  end
  ptr
end
