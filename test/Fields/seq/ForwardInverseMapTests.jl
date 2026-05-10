"""
Test the evaluation of the ForwardMap and InverseMap
To do: extend to 3D
"""

module ForwardInverseMapTests

using GridapGeosciences
using GridapGeosciences.Geometry
using Gridap
using Test

#### 2D test: four corners of the panel
pts_αβ = CUBE_HALF_EDGE.*[Point(-1.0,-1.0),Point(1.0,-1.0),Point(-1.0,1.0),Point(1.0,1.0)]

radius = 1

for panel in collect(1:NPANELS)
  fwd_maps = fill(ForwardMap(panel,radius),length(pts_αβ))
  inv_maps = fill(InverseMap(panel,radius),length(pts_αβ))

  pts_x = lazy_map(evaluate,fwd_maps,pts_αβ)

  pts_αβ_inv = lazy_map(evaluate,inv_maps,pts_x)

  @test pts_αβ ≈ pts_αβ
end




end # module
