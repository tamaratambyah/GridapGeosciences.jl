function reorder_cell_ids(model::DiscreteModel)
  n_cells_per_panel = Int(num_cells(model)/6)
  reindex = collect(1:n_cells_per_panel)
  return reindex
end

function reorder_cell_ids(model::Gridap.Adaptivity.AdaptedDiscreteModel{Dc}) where Dc
  println("no bang")
  n_cells_per_panel = Int(num_cells(model)/6)
  reindex = collect(1:n_cells_per_panel)

  reorder_cell_ids!(reindex,model.parent)

  return reindex
end

function reorder_cell_ids!(reindex,model::Gridap.Adaptivity.AdaptedDiscreteModel{Dc}) where Dc
  println("bang")

  n = Int(num_cells(model)/6)
  nc_fine = Tuple(fill(Int(sqrt(n)),Dc))
  f2c_cell_map,  = Gridap.Adaptivity._create_cartesian_f2c_maps(nc_fine, (2,2))
  c2f_cell_map = Adaptivity.get_o2n_faces_map(f2c_cell_map)
  p = collect(1: length(c2f_cell_map))

  reindex .= c2f_cell_map[p].data

  og_cell_map = copy(c2f_cell_map)

  reorder_cell_ids!(reindex,model.parent,c2f_cell_map,og_cell_map)

end

function reorder_cell_ids!(reindex,model::Gridap.Adaptivity.AdaptedDiscreteModel{Dc},
  prev_cell_map,og_cell_map) where Dc

  println("bang bang")

  n = Int(num_cells(model)/6)
  nc_fine = Tuple(fill(Int(sqrt(n)),Dc))
  _f2c_cell_map,  = Gridap.Adaptivity._create_cartesian_f2c_maps(nc_fine, (2,2))
  _c2f_cell_map = Adaptivity.get_o2n_faces_map(_f2c_cell_map)

  p = sortperm(_f2c_cell_map)
  o_reindex = prev_cell_map[p].data
  if length(p) == length(og_cell_map)
    reindex .= og_cell_map[p].data
  else
    reindex .= og_cell_map[o_reindex].data
  end

  reorder_cell_ids!(reindex,model.parent,_c2f_cell_map,og_cell_map)

end


function reorder_cell_ids!(reindex,model::DiscreteModel,c2f_cell_map,og_cell_map)
  println("coarest model")
end

function reorder_cell_ids!(reindex,model::DiscreteModel)
  println("coarest model")
end
