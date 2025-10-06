module PanelIdsTests
using Gridap
using GridapGeosciences
using Gridap.Adaptivity
using Test

################################################################################
## Test the panel ids from the BodyFittedTriangulation trian are the same as the
## panel ids from the model
################################################################################

function test_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  trian_panel_ids = get_panel_ids(Ω_panel)

  @test length(trian_panel_ids) == num_cells(panel_model)
  @test trian_panel_ids == get_panel_ids(panel_model)
end


s_panel_model = coarse_parametric_model()
test_panel_ids(s_panel_model)

s_panel_model = Gridap.Adaptivity.refine(s_panel_model)
test_panel_ids(s_panel_model)

s_panel_model = Gridap.Adaptivity.refine(s_panel_model)
test_panel_ids(s_panel_model)


end # module
