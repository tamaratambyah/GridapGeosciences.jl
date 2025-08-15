"""
The ambient model is no longer required.
This is because we now trick vtk into plotting on mapped nodes.
However, I will this here just inc ase
"""
function ambient_model(panel_model::ParametricDiscreteModel)
  panel_ids = get_panel_ids(panel_model)
  ambient_model(panel_model,panel_ids)
end

function ambient_model(panel_model,panel_ids)
  Dc = num_cell_dims(panel_model)
  Dp = num_point_dims(panel_model)

  @check Dc == Dp

  panel_grid = get_grid(panel_model)
  panel_topo = get_grid_topology(panel_model)
  ref_points = get_cell_ref_coordinates(panel_grid)

  ## map the panels to the ambient domain by applying the panel 1 forward map, and then rotation
  cmaps = get_cell_map(panel_grid)
  panel_to_ambient_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
  ambient_maps = lazy_map(∘,panel_to_ambient_map,cmaps)


  ## evaluate the ambient_nodes by applying the ambient_cmap to reference points
  ambient_cell_coords = lazy_map(evaluate,ambient_maps,ref_points)./1
  ambient_nodes = get_nodes_from_coords(panel_grid,ambient_cell_coords)

  ## the ambient_grid has the bespoke panel_2_ambient_cmap
  ambient_grid = Geometry.UnstructuredGrid(ambient_nodes,get_cell_node_ids(panel_grid),get_reffes(panel_grid),get_cell_type(panel_grid),OrientationStyle(panel_grid),
                      nothing,ambient_maps)
  ambient_topo = UnstructuredGridTopology(ambient_nodes,get_cell_node_ids(panel_grid),get_cell_type(panel_topo),get_polytopes(panel_topo),OrientationStyle(panel_topo))
  ambient_labels = FaceLabeling(ambient_topo)

  ambient_model = UnstructuredDiscreteModel(ambient_grid,ambient_topo,ambient_labels)

  return ambient_model
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
