panel_model = coarse_parametric_model()
panel_ids = get_panel_ids(panel_model)
cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
writevtk(Triangulation(panel_model),dir*"/ambient_grid_ref_lvl0",append=false,geo_map=cell_geo_map)

@check num_point_dims(panel_model) == num_cell_dims(panel_model) == 2

for n in [1,2]
  amodel = Gridap.Adaptivity.refine(panel_model)
  r_panel_ids = get_panel_ids(amodel)
  cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), r_panel_ids)
  writevtk(Triangulation(amodel),dir*"/ambient_grid_ref_lvl$(n)",append=false,geo_map=cell_geo_map)

  @check num_point_dims(amodel) == num_cell_dims(amodel) == 2

  @check is_child(amodel,panel_model)


  panel_model = amodel
end


nlevs = 4
panel_models = get_refined_models(nlevs)

for model in panel_models
  lvl = nref(nc(model))
  panel_ids = get_panel_ids(model)
  cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
  writevtk(Triangulation(model),dir*"/ambient_grid_ref_lvl$(lvl)",append=false,geo_map=cell_geo_map)
end




############# multigrid
using GridapSolvers
function ModelHierarchy(models::Vector{<:DiscreteModel})
  nlevs = length(models)
  @check all(map(i -> isa(models[i],AdaptedDiscreteModel),1:nlevs-1)) "Hierarchy models are not adapted models."
  for lev in 1:nlevs-1
    @check Adaptivity.is_child(models[lev],models[lev+1]) "Incorrect hierarchy of models."
  end

  level_parts = fill(DebugArray([1]),nlevs)
  meshes = Vector{ModelHierarchyLevel}(undef,nlevs)
  for lev in 1:nlevs-1
    glue  = Gridap.Adaptivity.get_adaptivity_glue(models[lev])
    model = models[lev]
    meshes[lev] = ModelHierarchyLevel(lev,model,glue,nothing,nothing)
  end
  meshes[nlevs] = ModelHierarchyLevel(nlevs,models[nlevs],nothing,nothing,nothing)
  return HierarchicalArray(meshes,level_parts)
end

using MPI, PartitionedArrays, GridapDistributed
import GridapSolvers.MultilevelTools: ModelHierarchyLevel, HierarchicalArray
mh = ModelHierarchy(panel_models)
th = Triangulation(mh)
