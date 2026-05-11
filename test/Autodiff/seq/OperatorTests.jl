"""
In this module, test that the sgrad, slap produce the same result as the
ambient_sgrad and ambient_slap
"""

module OperatorTests

using GridapGeosciences
using GridapGeosciences.Geometry
using Gridap
using Test

function test_operators(α,m)
  x = m(α)

  sgrad_ambient = ambient_sgrad(ambient_fX)(m)(x)
  sgrad_panel = sgrad(panel_fX)(m)(α)
  dif = sgrad_ambient .- sgrad_panel
  @test norm(dif) < 1e-12

  slap_ambient = ambient_surflap(ambient_fX)(m)(x)
  slap_panel = surflap(panel_fX)(m)(α)
  dif = slap_ambient .- slap_panel
  @test dif < 1e-12

end


function ambient_fX(xyz)
  x,y,z = xyz
  x^2 + y^2*z
end


function panel_fX(forward_map)
  function _f(αβ)
    x = forward_map(αβ)
    ambient_fX(x)
  end
end

panel_id = 1
radius = 1.0
thickness = 0.5

################################################################################
########## 2D
################################################################################
αβ = CUBE_HALF_EDGE*Point(-1.0,-1.0)
m = ForwardMap(panel_id,radius)
test_operators(αβ,m)

################################################################################
########## 3D
################################################################################
γαβ = Point(0.5,-CUBE_HALF_EDGE,-CUBE_HALF_EDGE)
m = ForwardMap(panel_id,radius,thickness)
test_operators(γαβ,m)


end # module
