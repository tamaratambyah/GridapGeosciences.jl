"""
In this module, we test the panel ids from the BodyFittedTriangulation trian are
equivalent to the panel ids from the model.
We also test the length of the panel ids is equvialent to the number of cells.
Complete this test for
  1. serial parametric model
  2. serial ambient model
"""

module PanelIdsTests
using Gridap
using GridapGeosciences
using Gridap.Adaptivity
using Test

function test_panel_ids(model)
  Ω_panel = Triangulation(model)
  trian_panel_ids = get_panel_ids(Ω_panel)

  @test length(trian_panel_ids) == num_cells(model)
  @test trian_panel_ids == get_panel_ids(model)
end

################################################################################
########## Parametric model
################################################################################
radius = 1.0
s_panel_model = coarse_parametric_model(radius)
test_panel_ids(s_panel_model)

s_panel_model = Gridap.Adaptivity.refine(s_panel_model)
test_panel_ids(s_panel_model)

s_panel_model = Gridap.Adaptivity.refine(s_panel_model)
test_panel_ids(s_panel_model)

################################################################################
########## Ambient model
################################################################################

s_ambient_model = CubedSphereAmbientDiscreteModel(radius;num_initial_uniform_refinements=0)
test_panel_ids(s_ambient_model)

s_ambient_model = Gridap.Adaptivity.refine(s_ambient_model)
test_panel_ids(s_ambient_model)

s_ambient_model = Gridap.Adaptivity.refine(s_ambient_model)
test_panel_ids(s_ambient_model)


@test true
end # module
