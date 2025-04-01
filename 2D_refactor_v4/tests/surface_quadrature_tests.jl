using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays


include("../src/initialise.jl")


function area_test(manifold_model,order::Int)
  computed_area = get_surface_area(manifold_model,order)
  _true_area = true_area(get_manifold_name(manifold_model))
  area_diff = abs(computed_area - _true_area)

  println("\nTrue area: $(_true_area)\nComputed area: $(computed_area)\nDifference: $(area_diff)")

  @test area_diff < 1e-8
end

## test the area of cube  -- apply 1 level of refinement
_manifold_model = ManifoldDiscreteModel(cube_model_3D,cube)
manifold_model = Adaptivity.refine(_manifold_model)

area_test(manifold_model,2)


## test the area of cubed sphere -- apply 1 level of refinement
_manifold_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
manifold_model = Adaptivity.refine(_manifold_model)
area_test(manifold_model,2)
area_test(manifold_model,4)
area_test(manifold_model,6)
area_test(manifold_model,8)




# quad = dΩg.s_quad.quad
# b = change_domain(f_cf,quad.trian,quad.data_domain_style)
# x = get_cell_points(quad)
# bx = b(x)

# g = m.metric
# _g = change_domain(g,quad.trian,quad.data_domain_style)
# x = get_cell_points(quad)
# gx = _g(x)
# gx_meas = lazy_map(MetricMeasure(),gx)

# y = bx[4].*gx_meas[4]
# z = lazy_map(Broadcasting(*), bx, gx_meas)
# println(y)
# println(z[4])
# println(y-z[4])
