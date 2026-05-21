using Gridap.Geometry
using Gridap.Arrays

function SingleAmbientDiscreteModel(panel_model::CartesianDiscreteModel,radius)
  model = UnstructuredDiscreteModel(panel_model)
  SingleAmbientDiscreteModel(model,radius)
end


function SingleAmbientDiscreteModel(panel_model::UnstructuredDiscreteModel,radius)

  panel_grid = get_grid(panel_model)
  panel_topo = get_grid_topology(panel_model)
  labels = get_face_labeling(panel_model)

  ## map: reffe -> alpha,beta
  cmap = get_cell_map(panel_grid)

  ## map: alpha,beta -> manifold
  fwd_map =  lazy_map(p -> ForwardMap(1,radius), collect(1:num_cells(panel_model)))

  ## map: reffe -> manifold
  geo_cmap = lazy_map(∘,fwd_map,cmap)

  ref_pts = get_cell_ref_coordinates(panel_grid)
  ambient_cell_coords = lazy_map(evaluate,geo_cmap,ref_pts)

  eT = eltype(testitem(ambient_cell_coords))
  panel_nodes = Gridap.Geometry.num_nodes(panel_grid)
  ambient_nodes = fill(zero(eT),panel_nodes)
  # Now create proper nodes from the cell-wise array of coords
  # We can do this for the ambient model because the nodes are defined uniquely
  # in ambient space
  cell_to_nodes = get_cell_node_ids(panel_grid)
  Gridap.FESpaces._free_and_dirichlet_values_fill!(
    ambient_nodes, eT[],
    array_cache(ambient_cell_coords),
    array_cache(cell_to_nodes),
    ambient_cell_coords,
    cell_to_nodes,
    eachindex(cell_to_nodes)
  )

  ## the ambient_grid has the bespoke panel_2_ambient_cmap
  ambient_grid = Gridap.Geometry.UnstructuredGrid(ambient_nodes,get_cell_node_ids(panel_grid),get_reffes(panel_grid),get_cell_type(panel_grid),OrientationStyle(panel_grid),
                      nothing,(geo_cmap))
  ambient_topo = UnstructuredGridTopology(ambient_nodes,get_cell_node_ids(panel_grid),get_cell_type(panel_topo),get_polytopes(panel_topo),OrientationStyle(panel_topo))
  ambient_labels = FaceLabeling(ambient_topo)

  UnstructuredDiscreteModel(ambient_grid,ambient_topo,labels)

end

function SingleParametricDiscreteModel(panel_model::CartesianDiscreteModel)
  SingleParametricDiscreteModel(UnstructuredDiscreteModel(panel_model))
end

function SingleParametricDiscreteModel(panel_model::UnstructuredDiscreteModel)

  panel_grid = get_grid(panel_model)
  panel_topo = get_grid_topology(panel_model)
  labels = get_face_labeling(panel_model)

  ## map: reffe -> alpha,beta
  cmap = get_cell_map(panel_grid)

  ref_pts = get_cell_ref_coordinates(panel_grid)
  ambient_cell_coords = lazy_map(evaluate,cmap,ref_pts)

  eT = eltype(testitem(ambient_cell_coords))
  panel_nodes = Gridap.Geometry.num_nodes(panel_grid)
  ambient_nodes = fill(zero(eT),panel_nodes)
  # Now create proper nodes from the cell-wise array of coords
  # We can do this for the ambient model because the nodes are defined uniquely
  # in ambient space
  cell_to_nodes = get_cell_node_ids(panel_grid)
  Gridap.FESpaces._free_and_dirichlet_values_fill!(
    ambient_nodes, eT[],
    array_cache(ambient_cell_coords),
    array_cache(cell_to_nodes),
    ambient_cell_coords,
    cell_to_nodes,
    eachindex(cell_to_nodes)
  )

  ## the ambient_grid has the bespoke panel_2_ambient_cmap
  ambient_grid = Gridap.Geometry.UnstructuredGrid(ambient_nodes,get_cell_node_ids(panel_grid),get_reffes(panel_grid),get_cell_type(panel_grid),OrientationStyle(panel_grid),
                      nothing,(cmap))
  ambient_topo = UnstructuredGridTopology(ambient_nodes,get_cell_node_ids(panel_grid),get_cell_type(panel_topo),get_polytopes(panel_topo),OrientationStyle(panel_topo))

  UnstructuredDiscreteModel(ambient_grid,ambient_topo,labels)

end
