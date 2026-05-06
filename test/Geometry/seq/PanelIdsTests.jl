"""
In this module, we test the panel ids from the BodyFittedTriangulation trian are
equivalent to the panel ids from the model.
We also test the length of the panel ids is equvialent to the number of cells.
"""

module PanelIdsTests
using Gridap
using GridapGeosciences
using Gridap.Adaptivity
using Test

function test_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  trian_panel_ids = get_panel_ids(Ω_panel)

  @test length(trian_panel_ids) == num_cells(panel_model)
  @test trian_panel_ids == get_panel_ids(panel_model)
end

radius = 1.0
s_panel_model = coarse_parametric_model(radius)
test_panel_ids(s_panel_model)

s_panel_model = Gridap.Adaptivity.refine(s_panel_model)
test_panel_ids(s_panel_model)

s_panel_model = Gridap.Adaptivity.refine(s_panel_model)
test_panel_ids(s_panel_model)

@test true
end # module
