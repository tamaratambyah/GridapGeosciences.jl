function parametric_model(cube_model)
  Dc = num_cell_dims(cube_model)
  Dp = num_point_dims(cube_model)

  @check Dp == Dc+1

  panel_ids = get_panel_ids(cube_model)
  cube_grid = get_grid(cube_model)
  cube_topo = get_grid_topology(cube_model)
  cube_nodes = get_node_coordinates(cube_grid)

  ## map the cube to the parametric domain
  cube_cmaps = get_cell_map(cube_grid)
  cube2panel_map = lazy_map(p->MatMultField( A_cube2panel[p] ), panel_ids)
  panel_cmaps = lazy_map(∘,cube2panel_map,cube_cmaps)


  ## create the panel model
  ## to correctly trigger Dc=2,Dp=2, need to have 2D nodes.
  ## the panel_nodes are just junk 2D nodes, never used by Gridap
  ## the panel_grid has the bespoke panel_cmap
  ## the panel_topo is the same as the cube, but with the 2D nodes so that Dp=2
  ## the panel_labels are from the panel_topo
  panel_nodes = map(x->Point(x[2],x[3]),cube_nodes) # these are just junk nodes, never used
  panel_grid = Geometry.UnstructuredGrid(panel_nodes,get_cell_node_ids(cube_grid),get_reffes(cube_grid),get_cell_type(cube_grid),OrientationStyle(cube_grid),
                      nothing,panel_cmaps)
  panel_topo = UnstructuredGridTopology(panel_nodes,get_cell_node_ids(cube_grid),get_cell_type(cube_topo),get_polytopes(cube_topo),OrientationStyle(cube_topo))
  panel_labels = FaceLabeling(panel_topo)
  panel_model = UnstructuredDiscreteModel(panel_grid,panel_topo,panel_labels)

  return panel_model,panel_ids
end
