function generate_ptr(n)
  nvertices = 4
  ptr  = Vector{Int}(undef,n+1)
  ptr[1]=1
  for i=1:n
    ptr[i+1]=ptr[i]+nvertices
  end
  ptr
end

function _CCAM_panel_wise_node_ids(npanels)
  ## CCAM panel ordering
  data = [ 1,2,3,4, 3,4,5,6,  2,7,4,6, 8,5,7,6, 1,8,2,7, 1,3,8,5 ]
  ptr = generate_ptr(npanels)
  Table(data,ptr)
end 

function _CCAM_cube_nodes_3d(a::Real)
   a.* [
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


function coarse_cube_surface_3D(a::Real,npanels::Int)

  nodes_3d = _CCAM_cube_nodes_3d(a)
  cell_node_ids = _CCAM_panel_wise_node_ids(npanels)

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



function coarse_cube_model(a::Real,npanels::Int)
  cube_grid,topo,labels, = coarse_cube_surface_3D(a,npanels)
  model = UnstructuredDiscreteModel(cube_grid,topo,labels)
  return model
end
