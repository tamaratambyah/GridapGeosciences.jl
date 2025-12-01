

omodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=1)
panel_model = omodel.parametric_dmodel
panel_ids = get_panel_ids(panel_model)

mask = map(panel_ids) do p
  return p .== pid
end

gids = get_cell_gids(panel_model)
map(partition(gids)) do gid
  own_to_local(gid)
  local_to_ghost(gid)
end

models = map(local_views(panel_model),mask) do panel_model, mask
  println(findall(mask))
  model = DiscreteModelPortion(panel_model,mask)
  _grid = get_grid(model)

  ## construct the new grid by hand, so that specific cmaps are given
  nodes = collect1d(Geometry.get_node_coordinates(_grid)) # these are just junk nodes, never used
  ctype = collect1d(get_cell_type(_grid))
  cmaps = collect1d(get_cell_map(_grid))
  grid = Gridap.Geometry.UnstructuredGrid(nodes,get_cell_node_ids(_grid),
            get_reffes(_grid),ctype,OrientationStyle(_grid),
                      nothing,cmaps)

  topo = UnstructuredGridTopology(get_grid_topology(model))
  labels = FaceLabeling(topo)

  UnstructuredDiscreteModel(grid,topo,labels)
end

metadata = nothing
single_panel_model = GridapDistributed.GenericDistributedDiscreteModel(models,get_cell_gids(panel_model);metadata)
