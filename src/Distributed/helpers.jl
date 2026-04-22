#### The distributed panel ids are extracted from the serial. This includes both
#### owned+ghost panel_ids.
#### The owned panel_ids are extracted by determining the owned cell
function distributed_panel_ids(dmodel,spanel_ids::AbstractArray{Int})
  gids = get_cell_gids(dmodel)

  dpanel_ids = map(partition(gids)) do ids
    lid_to_gid = local_to_global(ids)
    return spanel_ids[lid_to_gid]
  end

  owned_panel_ids = map(dpanel_ids,partition(gids)) do panel_ids, cids
    owned_cells = own_to_local(cids)
    return panel_ids[owned_cells]
  end

  return dpanel_ids, owned_panel_ids
end


### This is Jordi's function from GridapDistributed/test/AdaptivityTests.jl
function DistributedAdaptivityGlue(serial_glue,parent,child)
  glue = map(partition(get_cell_gids(parent)),partition(get_cell_gids(child))) do parent_gids, child_gids
    old_g2l = global_to_local(parent_gids)
    old_l2g = local_to_global(parent_gids)
    new_l2g = local_to_global(child_gids)

    n2o_cell_map  = lazy_map(Reindex(old_g2l),serial_glue.n2o_faces_map[3][new_l2g])
    n2o_faces_map = [Int64[],Int64[],collect(n2o_cell_map)]
    n2o_cell_to_child_id = serial_glue.n2o_cell_to_child_id[new_l2g]
    rrules = serial_glue.refinement_rules[old_l2g]
    Gridap.Adaptivity.AdaptivityGlue(n2o_faces_map,n2o_cell_to_child_id,rrules)
  end
  return glue
end




### return distributed version of serial model after n_ref_lvls of refinement
function get_distributed_panel_model(ranks,nprocs,n_ref_lvls::Int,radius)

  # get refined models in serial
  s_panel_model = coarse_parametric_model(radius)
  for n in 1:n_ref_lvls
    s_panel_model = Gridap.Adaptivity.refine(s_panel_model)
  end

  spanel_ids = get_panel_ids(s_panel_model)

  part_to_cells = [PartitionedArrays.local_range(rank,nprocs,num_cells(s_panel_model)) for rank in 1:nprocs]

  # get the partition
  fine_cell_to_part = zeros(Int32,num_cells(s_panel_model))
  for (rank, cells) in enumerate(part_to_cells)
    fine_cell_to_part[cells] .= rank
  end

  # distribute the model
  dmodel = UnstructuredDiscreteModel(DiscreteModel(ranks,Adaptivity.get_model(s_panel_model),fine_cell_to_part))
  dpanel_ids, = distributed_panel_ids(dmodel,spanel_ids)
  panel_model = DistributedParametricDiscreteModel(dmodel,dpanel_ids,radius)
 return panel_model
end
