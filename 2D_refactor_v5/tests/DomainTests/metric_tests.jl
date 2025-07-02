"""
Test the Jacobian transformation of each panel is equivalent to g=J^tJ
The Jacobian is of the mapping κ: α → X
"""
using DrWatson
using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Fields
using Test
include("../../src/initialise.jl")

function lazy_test(d)
  for i in eachindex(d)
    @test sum( norm.(d[i]) .< 1e-14 ) == length(d[i])
  end
  println("Yay!")
end

manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(manifold_model)
panel_ids = get_panel_ids(manifold_model)
Ω = Triangulation(manifold_model)
pts = get_cell_points(Ω)
m = Metric(cubedsphere,Ω)

"""
check panel wise jacobians are the rotation of the panel 1 jacobian
"""
Jp(p) = x -> r1p_3D[p]⋅Jpanel1(x)

Jp(2)(Point(1,1)) == Jpanel2(Point(1,1))
Jp(3)(Point(1,1)) == Jpanel3(Point(1,1))
Jp(4)(Point(1,1)) == Jpanel4(Point(1,1))
Jp(5)(Point(1,1)) == Jpanel5(Point(1,1))
Jp(6)(Point(1,1)) == Jpanel6(Point(1,1))


"""
Compute the metric, inverse and sqrt(meas) using Jacobians
"""
_Jcf = map(p->GenericField(Jp(p)),panel_ids)
Jcf = CellData.GenericCellField(_Jcf,Ω,PhysicalDomain())
Gcf = (Operation(transpose)(Jcf)) ⋅ Jcf
invG = Operation(inv)(Gcf)
measG = Operation((meas))(Gcf)
sqrtmeasG = Operation(sqrt)(measG)


lazy_test((Gcf - m.metric)(pts))
lazy_test( (invG - m.inv_metric)(pts) )
lazy_test( (sqrtmeasG-m.sq_meas)(pts) )

"""
compute sqrt g analytically
"""
function sqrtgfunc(αβ)
  α,β = αβ
  RADIUS^2/( (1 + (tan(α))^2 + (tan(β))^2 )^(3/2)  *(cos(α))^2*(cos(β))^2 )
end
sqrtg = CellField(sqrtgfunc,Ω)
lazy_test( (sqrtmeasG - sqrtg)(pts) )
lazy_test( (m.sq_meas - sqrtg)(pts) )
