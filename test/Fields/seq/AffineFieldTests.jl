
"""
Test the evaluation of the AffineField that performs the operation Ax + b
"""

module AffineFieldTests

using Gridap
using GridapGeosciences
using Test

using GridapGeosciences.Geometry
using GridapGeosciences.Fields

panel_ids = collect(1:NPANELS)
A_panel2cube_transpose = map(x->transpose(x),A_panel2cube)

panel2cube_map = lazy_map(p-> MyAffineField(A_panel2cube[p],b_panel2cube[p]), panel_ids)
∇panel2cube_map = lazy_map(gradient,panel2cube_map)

# Test evaluation on a single point
x = Point(-CUBE_HALF_EDGE,-CUBE_HALF_EDGE)
evaluate(panel2cube_map,x)
@test evaluate(∇panel2cube_map,x) == A_panel2cube_transpose

## Repeat for array of points
xt = fill(x,length(panel_ids))
lazy_map(evaluate,panel2cube_map,xt)
@test lazy_map(evaluate,∇panel2cube_map,xt) == A_panel2cube_transpose

@test true


end # module
