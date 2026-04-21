"""
In this module, we test the refinement of the serial model by checking
1. the number of cell dims and point dims = 2
2. the refined model is a child of the parent
"""

module RefinementTests

using Gridap
using GridapGeosciences
using Gridap.Helpers,  Gridap.Adaptivity
using Test

### Check the Dc, Dp of the coarse model
radius = 1.0
panel_model = coarse_parametric_model(radius)
panel_ids = get_panel_ids(panel_model)
# writevtk(Triangulation(panel_model),"ambient_grid_ref_lvl0",append=false,geo_map=geo_map_func(panel_ids))

@test num_point_dims(panel_model) == num_cell_dims(panel_model) == 2

### Apply refinement, check the list of refined models
n_ref_lvls = 4
panel_models = get_refined_models(n_ref_lvls, radius)
for lev in 1:n_ref_lvls-1
  @test num_point_dims(panel_models[lev]) == num_cell_dims(panel_models[lev]) == 2
  @test is_child(panel_models[lev],panel_models[lev+1])
end

@test true
end
