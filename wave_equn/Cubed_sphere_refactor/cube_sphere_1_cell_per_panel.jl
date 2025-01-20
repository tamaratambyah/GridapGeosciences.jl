"""
cube_sphere_1_cell_per_panel -- builds the flat cube geometry, with correct
periodicity, with 1 cell per panel based on the CCAM ordering of panels. This
geometry contains 8 nodes. The following table information is required:

panel/cell no | cell_node_ids
      1       |   [1 2 3 4]
      2       |   [3 4 5 6]
      3       |   [4 2 6 7]
      4       |   [6 7 5 8]
      5       |   [7 2 8 1]
      6       |   [8 1 5 3]

  node no     | Point in 3D
      1       |   (1,-1,-1)
      2       |   (1,1,-1)
      3       |   (1,-1,1)
      4       |   (1,1,1)
      5       |   (-1,-1,1)
      6       |   (-1,1,1)
      7       |   (-1,1,-1)
      8       |   (-1,-1,-1)

"""
function cube_sphere_1_cell_per_panel()
  npanels = 6

  nodes = [
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
  data = [ 1,2,3,4, 3,4,5,6, 4,2,6,7, 6,7,5,8, 7,2,8,1, 8,1,5,3  ]
  ptr = generate_ptr(npanels)
  cell_node_ids = Table(data,ptr)

  polytopes = fill(QUAD,npanels)
  cell_type = fill(1,npanels)
  reffes = LagrangianRefFE(Float64,QUAD,1)
  cell_reffes=[reffes]

  topo = UnstructuredGridTopology(nodes,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
  face_labels = FaceLabeling(topo)

  cube_grid = Gridap.Geometry.UnstructuredGrid(nodes,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

  return cube_grid,topo,face_labels
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
