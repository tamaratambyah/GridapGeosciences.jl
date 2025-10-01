module RefinementTests

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


end
