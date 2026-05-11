"""
In this module, test that the sgrad, slap produce the same result as the
ambient_sgrad and ambient_slap
"""

module OperatorTests

using GridapGeosciences
using GridapGeosciences.Geometry
using Gridap
using Test

function ambient_fX(xyz)
  x,y,z = xyz
  x*y*z
end


function panel_fX(forward_map)
  function _f(αβ)
    x = forward_map(αβ)
    ambient_fX(x)
  end
end

panel_id = 1
radius = 1


αβ = CUBE_HALF_EDGE*Point(-1.0,-1.0)

m = ForwardMap(panel_id,radius)
minv = InverseMap(m)

x = m(αβ)

sgrad_ambient = ambient_sgrad(ambient_fX)(m)(x)
sgrad_panel = sgrad(panel_fX)(m)(αβ)
dif = sgrad_ambient .- sgrad_panel
@test norm(dif) < 1e-12

slap_ambient = ambient_surflap(ambient_fX)(m)(x)
slap_panel = surflap(panel_fX)(m)(αβ)
dif = slap_ambient .- slap_ambient
@test dif < 1e-12




end # module
