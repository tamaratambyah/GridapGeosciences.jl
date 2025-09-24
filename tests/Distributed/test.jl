using DrWatson
using Gridap
using GridapDistributed
using PartitionedArrays
using MPIPreferences

using Gridap.Geometry
using Gridap.Adaptivity
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

function distributed_pids(model,spids)
  gids = get_cell_gids(model)
  pids = map(partition(gids)) do ids
    lid_to_gid = local_to_global(ids)
    return spids[lid_to_gid]
  end
  return pids
end

################################################################################
using GridapGeosciences
s_model_coarse = coarse_parametric_model()
s_model_ref = refine(s_model_coarse)
s_model_ref_ref = refine(s_model_ref)


# serial_panel_ids = get_panel_ids(s_model_coarse)
# serial_cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), serial_panel_ids)
# writevtk(Triangulation(s_model_coarse),dir*"/ambient_model",geo_map=serial_cell_geo_map)

s_model_coarse = Geometry.UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(2,3)))
s_model_ref = refine(s_model_coarse)
s_model_ref_ref = refine(s_model_ref)

################################################################################
##### Distributed models
################################################################################
nprocs = (2,2)
ranks  = with_debug() do distribute
  distribute(LinearIndices((prod(nprocs),)))
end

## Level 1 of refinement
#parent_cell_to_part = sort(rand(1:4,num_cells(s_model_coarse))) #[1,1,2,2,3,4]

part_to_cells = [PartitionedArrays.local_range(rank,prod(nprocs),num_cells(s_model_coarse)) for rank in 1:prod(nprocs)]
coarse_cell_to_part = zeros(Int32,num_cells(s_model_coarse))
for (rank, cells) in enumerate(part_to_cells)
  coarse_cell_to_part[cells] .= rank
end

spids = [get_panel_ids(s_model_ref_ref),get_panel_ids(s_model_ref),get_panel_ids(s_model_coarse)]

models = [s_model_ref_ref.model,s_model_ref.model,s_model_coarse]
glues = [s_model_ref_ref.glue,s_model_ref.glue]
cell_to_part = Vector{Any}(undef,length(models))
cell_to_part[end] = coarse_cell_to_part
for level in length(models)-1:-1:1
  n2o_cells = glues[level].n2o_faces_map[3]
  cell_to_part[level] = cell_to_part[level+1][n2o_cells]
end

dmodels = Vector{Any}(undef,length(models))
pids = Vector{Any}(undef,length(models))
dmodels[end] = DiscreteModel(ranks,models[end],cell_to_part[end])
pids[end] = distributed_pids(dmodels[end],spids[end])
for level in length(models)-1:-1:1
  child = Gridap.Geometry.UnstructuredDiscreteModel(DiscreteModel(ranks,models[level],cell_to_part[level]))
  parent = dmodels[level+1]
  glue = DistributedAdaptivityGlue(glues[level],parent,child)
  dmodels[level] = GridapDistributed.DistributedAdaptedDiscreteModel(child,parent,glue)
  pids[level] = distributed_pids(child,spids[level])
end

for (level,(model,cparts,panel_ids)) in enumerate(zip(models,cell_to_part,spids))

  #include("vtk.jl")
  geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
  writevtk(Triangulation(model),dir*"/serial_model_ref_$(level)";append=false,celldata=["part" => cparts],geo_map=geo_map)

end

include("vtk.jl")
for (level,(_dmodel,panel_ids)) in enumerate(zip(dmodels,pids))
  dmodel = (level == 3) ? _dmodel : Gridap.Adaptivity.get_model(_dmodel)

  gids = get_cell_gids(dmodel)
  cell_geo_map = map(panel_ids,partition(gids)) do panel_ids, cids
    owned_cells = own_to_local(cids)
    println(panel_ids[owned_cells])
    return lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids[owned_cells])
  end

  writevtk(Triangulation(dmodel),dir*"/ambient_model_ref_$(level)"; append=false)#, geo_map=cell_geo_map)
end

child_cell_to_part  = lazy_map(Reindex(parent_cell_to_part),s_ref_glue.n2o_faces_map[3])
println(child_cell_to_part)
if GridapDistributed.i_am_in(ranks)
  parent = DiscreteModel(ranks,s_model_coarse,parent_cell_to_part)
  child  = DiscreteModel(ranks,s_model_ref,child_cell_to_part)
  coarse_adaptivity_glue = DistributedAdaptivityGlue(s_ref_glue,parent,child)
else
  println("not main ")
  parent = nothing; child  = nothing; coarse_adaptivity_glue = nothing
end

model_coarse = parent
model_ref = GridapDistributed.DistributedAdaptedDiscreteModel(child,parent,coarse_adaptivity_glue)
# Gridap.Adaptivity.is_child(model_ref,parent)

## Level 2 of refinement
ref_parent_cell_to_part = sort(rand(1:4,num_cells(s_model_ref)))
ref_child_cell_to_part  = lazy_map(Reindex(ref_parent_cell_to_part),s_ref_ref_glue.n2o_faces_map[3])
println(ref_child_cell_to_part)
if GridapDistributed.i_am_in(ranks)
  ref_parent = DiscreteModel(ranks,s_model_ref,ref_parent_cell_to_part) # DiscreteModel(ranks,serial_parent,parent_cell_to_part)
  ref_child  = DiscreteModel(ranks,s_model_ref_ref,ref_child_cell_to_part) #DiscreteModel(ranks,serial_child,child_cell_to_part)
  ref_adaptivity_glue = DistributedAdaptivityGlue(s_ref_ref_glue,ref_parent,ref_child)
else
  ref_parent = nothing; ref_child  = nothing; ref_adaptivity_glue = nothing
end

# model_ref = parent
model_ref_ref = GridapDistributed.DistributedAdaptedDiscreteModel(ref_child,ref_parent,ref_adaptivity_glue)
Gridap.Adaptivity.is_child(model_ref_ref,ref_parent)


# model_coarse = parent
# model_ref = ref_parent
# model_ref_ref = ref_child


########## Coarse model -- distributed form
### get the serial panel ids and map to distributed cells
serial_panel_ids = get_panel_ids(s_model_coarse)
dmodel = model_coarse

mgids = get_cell_gids(dmodel)

notcells, tcell_to_mcell = map(
      local_views(dmodel),local_views(dtrian),partition(mgids)) do model,trian,partition
      lid_to_owner = local_to_owner(partition)
      part = part_id(partition)
      glue = get_glue(trian,Val(2))
      @assert isa(glue,Gridap.Geometry.FaceToFaceGlue)
      tcell_to_mcell = glue.tface_to_mface
      notcells = count(tcell_to_mcell) do mcell
        lid_to_owner[mcell] == part
      end
      notcells, tcell_to_mcell
    end |> tuple_of_arrays

# Find the global range of owned dofs
first_gtcell = scan(+,notcells,type=:exclusive,init=one(eltype(notcells)))

mcell_to_gtcell = map(
  first_gtcell,tcell_to_mcell,partition(mgids)) do first_gtcell,tcell_to_mcell,partition
  mcell_to_gtcell = zeros(Int,local_length(partition))
  loc_to_owner = local_to_owner(partition)
  part = part_id(partition)
  gtcell = first_gtcell
  for mcell in tcell_to_mcell
    if loc_to_owner[mcell] == part
      mcell_to_gtcell[mcell] = serial_panel_ids[gtcell] ## extract the serial panel id
      gtcell += 1
    end
  end
  mcell_to_gtcell
end

cell_geo_map = map(mcell_to_gtcell) do pid
  pids = filter(x -> x > 0, pid)
  return lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), pids)
end

include("vtk.jl")
writevtk(Triangulation(dmodel),dir*"/ambient_model_ref0",geo_map=cell_geo_map)


########## Level 1 refined -- distributed form
serial_panel_ids = get_panel_ids(s_model_ref)
dmodel = model_ref
dtrian = Triangulation(dmodel)
mgids = get_cell_gids(dmodel)

notcells, tcell_to_mcell = map(
      local_views(dmodel),local_views(dtrian),partition(mgids)) do model,trian,partition
      lid_to_owner = local_to_owner(partition)
      part = part_id(partition)
      glue = get_glue(trian,Val(2))
      @assert isa(glue,Gridap.Geometry.FaceToFaceGlue)
      tcell_to_mcell = glue.tface_to_mface
      notcells = count(tcell_to_mcell) do mcell
        lid_to_owner[mcell] == part
      end
      notcells, tcell_to_mcell
    end |> tuple_of_arrays

# Find the global range of owned dofs
first_gtcell = scan(+,notcells,type=:exclusive,init=one(eltype(notcells)))



mcell_to_gtcell = map(
  first_gtcell,tcell_to_mcell,partition(mgids)) do first_gtcell,tcell_to_mcell,partition
  mcell_to_gtcell = zeros(Int,local_length(partition))
  loc_to_owner = local_to_owner(partition)
  part = part_id(partition)
  gtcell = first_gtcell
  for mcell in tcell_to_mcell
    if loc_to_owner[mcell] == part
      mcell_to_gtcell[mcell] = serial_panel_ids[gtcell] ## extract the serial panel id
      gtcell += 1
    end
  end
  mcell_to_gtcell
end


cell_geo_map = map(mcell_to_gtcell) do pid
  pids = filter(x -> x > 0, pid)
  return lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), pids)
end

include("vtk.jl")
writevtk(Triangulation(dmodel),dir*"/ambient_model_ref1",geo_map=cell_geo_map)


######### another level of refinement
serial_panel_ids = get_panel_ids(s_model_ref_ref)
dmodel = model_ref_ref
dtrian = Triangulation(dmodel)
mgids = get_cell_gids(dmodel)

notcells, tcell_to_mcell = map(
      local_views(dmodel),local_views(dtrian),partition(mgids)) do model,trian,partition
      lid_to_owner = local_to_owner(partition)
      part = part_id(partition)
      glue = get_glue(trian,Val(2))
      @assert isa(glue,Gridap.Geometry.FaceToFaceGlue)
      tcell_to_mcell = glue.tface_to_mface
      notcells = count(tcell_to_mcell) do mcell
        lid_to_owner[mcell] == part
      end
      notcells, tcell_to_mcell
    end |> tuple_of_arrays

# Find the global range of owned dofs
first_gtcell = scan(+,notcells,type=:exclusive,init=one(eltype(notcells)))

mcell_to_gtcell = map(
  first_gtcell,tcell_to_mcell,partition(mgids)) do first_gtcell,tcell_to_mcell,partition
  mcell_to_gtcell = zeros(Int,local_length(partition))
  loc_to_owner = local_to_owner(partition)
  part = part_id(partition)
  gtcell = first_gtcell
  for mcell in tcell_to_mcell
    if loc_to_owner[mcell] == part
      mcell_to_gtcell[mcell] = serial_panel_ids[gtcell] ## extract the serial panel id
      gtcell += 1
    end
  end
  mcell_to_gtcell
end


cell_geo_map = map(mcell_to_gtcell) do pid
  pids = filter(x -> x > 0, pid)
  return lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), pids)
end

writevtk(Triangulation(dmodel),dir*"/ambient_model_ref2",geo_map=cell_geo_map)
