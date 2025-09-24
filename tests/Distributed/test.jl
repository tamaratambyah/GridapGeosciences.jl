using DrWatson
using Gridap
using GridapDistributed
using PartitionedArrays
using MPIPreferences
MPIPreferences.use_jll_binary()

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

################################################################################
using GridapGeosciences
serial_parent = coarse_parametric_model()
serial_child = refine(serial_parent)
serial_rglue = Gridap.Adaptivity.get_adaptivity_glue(serial_child)

serial_panel_ids = get_panel_ids(serial_parent)
serial_cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), serial_panel_ids)
writevtk(Triangulation(serial_parent),dir*"/ambient_model",geo_map=serial_cell_geo_map)



nprocs = (1,2)
ranks  = with_debug() do distribute
  distribute(LinearIndices((prod(nprocs),)))
end

parent_cell_to_part = [1,1,1,2,2,2]
child_cell_to_part  = lazy_map(Reindex(parent_cell_to_part),serial_rglue.n2o_faces_map[3])
println(child_cell_to_part)
if GridapDistributed.i_am_in(ranks)
  parent = DiscreteModel(ranks,serial_parent,parent_cell_to_part)
  child  = DiscreteModel(ranks,serial_child,child_cell_to_part)
  coarse_adaptivity_glue = DistributedAdaptivityGlue(serial_rglue,parent,child)
else
  parent = nothing; child  = nothing; coarse_adaptivity_glue = nothing
end

model_ref = GridapDistributed.DistributedAdaptedDiscreteModel(child,parent,coarse_adaptivity_glue)
Gridap.Adaptivity.is_child(model_ref,parent)

#### find owned cells
cell_gids = get_cell_gids(parent)

map(partition(cell_gids)) do gid
  # local_to_owner(gid)
  # local_to_global(gid)
  part_id(gid)
end

##### Try to plot in vtk
panel_ids = map(local_views(parent)) do model
  println(typeof(model))
  return get_panel_ids(model)
end

cell_geo_map = map(panel_ids) do pid
  for p in pid
    println(p)
  end
  return lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), pid)
end

include("vtk.jl")
writevtk(Triangulation(parent),dir*"/ambient_model",geo_map=cell_geo_map)

# map(local_views(model_ref)) do model
#   println(typeof(model))
#   get_panel_ids(model)
# end

cellx = get_cell_coordinates(trian)
  println(cellx)
  cellx_mapped = lazy_map(evaluate,geo_map,cellx)







################################################################################
##### debug
################################################################################
model = CartesianDiscreteModel(ranks,nprocs,(0,1,0,1),(3,2))
model_ref = GridapDistributed.Adaptivity.refine(model)
model_ref_ref = GridapDistributed.Adaptivity.refine(model_ref)
model_ref_ref_ref = GridapDistributed.Adaptivity.refine(model_ref_ref)
# get_panel_ids!(panel_ids,model)

# panel_ids = map(local_views(model) ) do model
#   collect(1:num_cells(model))
# end
# get_panel_ids!(panel_ids,model)
# panel_ids

Dc = 2

panel_ids = map(local_views(model_ref_ref)) do model
  glue = Gridap.Adaptivity.get_adaptivity_glue(model)
  f2c_map = glue.n2o_faces_map[Dc+1]
  return f2c_map
end

map(GridapDistributed.Adaptivity.get_parent,local_views(model_ref_ref)) do local_parent
  println(typeof(local_parent))
end

parent = GridapDistributed.Adaptivity.get_parent(model_ref_ref)


map(local_views(parent)) do model
  glue = Gridap.Adaptivity.get_adaptivity_glue(model)
  f2c_map = glue.n2o_faces_map[Dc+1]
  panel_ids .= n2o[panel_ids]
  return panel_ids
end





function get_panel_ids(model::GridapDistributed.DistributedAdaptedDiscreteModel{Dc}) where Dc
  # panel_ids = copy(model.glue.n2o_faces_map[Dc+1])
  panel_ids = map(local_views(model_ref)) do model
    glue = Gridap.Adaptivity.get_adaptivity_glue(model)
    f2c_map = glue.n2o_faces_map[Dc+1]
    return f2c_map
  end

  get_panel_ids!(panel_ids, model.parent)
  return panel_ids
end


function get_panel_ids!(panel_ids, model::GridapDistributed.DistributedAdaptedDiscreteModel{Dc}) where Dc
  n2o = model.glue.n2o_faces_map[Dc+1]
  panel_ids .= n2o[panel_ids]
  get_panel_ids!(panel_ids, model.parent)
end



function get_panel_ids(model::GridapDistributed.DistributedDiscreteModel)
  panel_ids = map(local_views(model) ) do model
    collect(1:num_cells(model))
  end
  return panel_ids
end

function get_panel_ids!(panel_ids, model::GridapDistributed.DistributedDiscreteModel)
end


panel_model = coarse_parametric_model()
