using DrWatson
using Gridap
using GridapGeosciences
using Gridap.Helpers,  Gridap.Adaptivity
using Test

### Check the Dc, Dp of the coarse model
panel_model = coarse_parametric_model()
panel_ids = get_panel_ids(panel_model)
cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
writevtk(Triangulation(panel_model),dir*"/ambient_grid_ref_lvl0",append=false,geo_map=cell_geo_map)

@test num_point_dims(panel_model) == num_cell_dims(panel_model) == 2

### Apply refinement
n_ref_lvls = 4
panel_model = coarse_parametric_model()
for n in collect(1:n_ref_lvls)
  amodel = refine(panel_model)
  r_panel_ids = get_panel_ids(amodel)
  cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), r_panel_ids)
  writevtk(Triangulation(amodel),dir*"/ambient_grid_ref_lvl$(n)",append=false,geo_map=cell_geo_map)

  @test num_point_dims(amodel) == num_cell_dims(amodel) == 2
  @test is_child(amodel,panel_model)

  panel_model = amodel
end


### Check the list of refined models
panel_models = get_refined_models(n_ref_lvls)
for lev in 1:n_ref_lvls-1
  @test is_child(panel_models[lev],panel_models[lev+1])
end



############# multigrid
using GridapSolvers
using GridapDistributed, PartitionedArrays

model0 = coarse_parametric_model()
model1 = Adaptivity.refine(model0)
model2 = Adaptivity.refine(model1)
models = [model2,model1]
mh = GridapSolvers.MultilevelTools.ModelHierarchy(models)


th = Triangulation(mh)



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
