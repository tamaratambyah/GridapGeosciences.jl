### refine near the boundary of domain [0,1]^2
function boundary_refinement(dmodel::GridapDistributed.DistributedDiscreteModel)
  ref_coarse_flags=map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
    flags=zeros(Cint,length(indices))
    flags.=nothing_flag

    flags[1]=refine_flag
    flags[own_length(indices)]=refine_flag

    # To create some unbalance
    if (rank%2==0 && own_length(indices)>1)
        flags[div(own_length(indices),2)]=refine_flag
    end
    flags
  end
end

### refine in the middle of the domain [0,1]^2
function middle_refinement(dmodel::GridapDistributed.DistributedDiscreteModel)

  ref_coarse_flags=map(partition(get_cell_gids(dmodel.dmodel)), local_views(dmodel) ) do indices,lmodel
    flags=zeros(Cint,length(indices))
    flags.=nothing_flag

    cmap = get_cell_map(get_grid(lmodel))
    ref_points = get_cell_ref_coordinates(lmodel)
    coords = lazy_map(evaluate,cmap,ref_points)

    for (i,xy) in enumerate(coords)
      x = map(x->x[1],xy)
      y = map(x->x[2],xy)
      if any(x .> 0.4 ) && any(x .< 0.6 ) && any(y .> 0.4) && any(y .< 0.6)
          flags[i] = refine_flag
      end
    end
    flags
  end
  return ref_coarse_flags
end


function initial_unbalance(dmodel::GridapDistributed.DistributedDiscreteModel)

  ref_coarse_flags=map(partition(get_cell_gids(dmodel.dmodel)), local_views(dmodel) ) do indices,lmodel
    flags=zeros(Cint,length(indices))
    flags.=nothing_flag

    cmap = get_cell_map(get_grid(lmodel))
    ref_points = get_cell_ref_coordinates(lmodel)
    coords = lazy_map(evaluate,cmap,ref_points)

    for (i,xy) in enumerate(coords)
      x = map(x->x[1],xy)
      y = map(x->x[2],xy)
      if any(x .> 0.4 ) && any(x .< 0.6 ) && any(y .> 0.4) && any(y .< 0.6)
          flags[i] = refine_flag
      end
      if any(x .< 0.2 )  && any(y .< 0.2)
        flags[i] = refine_flag
      end
      if any(x .> 0.8 )  && any(y .> 0.8)
        flags[i] = refine_flag
      end

    end
    flags
  end
  return ref_coarse_flags
end


function refine_all(model::GridapDistributed.DistributedDiscreteModel)
  cell_partition=get_cell_gids(model)
  ref_coarse_flags=map(partition(cell_partition)) do indices
    flags=zeros(Cint,length(indices))
    flags.=refine_flag
  end
  fmodel, = Gridap.Adaptivity.adapt(model,ref_coarse_flags)
  return fmodel
end
