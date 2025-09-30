
function get_distributed_refined_models(ranks,nprocs,n_ref_lvls::Int,coarse_s_model=true)

  # get refined models in serial
  s_models  = get_refined_models(n_ref_lvls,coarse_s_model)
  spanel_ids = map(m->get_panel_ids(m),s_models)
  s_model_coarse = s_models[end]

  # extract the models and glues in arrays
  models = map(m->get_model(m),s_models[1:end-1])
  glues = map(m->get_adaptivity_glue(m),s_models[1:end-1])
  if coarse_s_model
    push!(models,s_model_coarse)
  else
    push!(models,get_model(s_model_coarse))
  end

  # partition the processors
  part_to_cells = [PartitionedArrays.local_range(rank,nprocs,num_cells(s_model_coarse)) for rank in 1:nprocs]

  # store the partition for each model
  cell_to_part = Vector{Vector{Int32}}(undef,length(models))

  # get the coarse partition
  coarse_cell_to_part = zeros(Int32,num_cells(s_model_coarse))
  for (rank, cells) in enumerate(part_to_cells)
    coarse_cell_to_part[cells] .= rank
  end
  cell_to_part[end] = coarse_cell_to_part

  # get the refine partition based on glue
  for level in length(models)-1:-1:1
    n2o_cells = glues[level].n2o_faces_map[3]
    cell_to_part[level] = cell_to_part[level+1][n2o_cells]
  end

  # construct array of distributed models, distributed panel ids (all + owned only)
  dmodels = Vector{DistributedParametricDiscreteModel}(undef,length(models))
  dpanel_ids = Vector{AbstractArray{Vector{Int}}}(undef,length(models))
  owned_panel_ids = Vector{AbstractArray{Vector{Int}}}(undef,length(models))

  # the coarsest model is the last in the list
  coarse_dmodel = DiscreteModel(ranks,models[end],cell_to_part[end])
  dpanel_ids[end],owned_panel_ids[end] = distributed_panel_ids(coarse_dmodel,spanel_ids[end])
  dmodels[end] = DistributedParametricDiscreteModel(coarse_dmodel,dpanel_ids[end])

  # loop backwards through refinement levels
  for level in length(models)-1:-1:1
    child = DiscreteModel(ranks,models[level],cell_to_part[level])
    parent = dmodels[level+1]
    glue = DistributedAdaptivityGlue(glues[level],parent,child)
    dpanel_ids[level],owned_panel_ids[level] = distributed_panel_ids(child,spanel_ids[level])

    dmodels[level] = DistributedAdaptedParametricDiscreteModel(child,parent,glue,dpanel_ids[level])
  end

  return dmodels, dpanel_ids, owned_panel_ids

end
