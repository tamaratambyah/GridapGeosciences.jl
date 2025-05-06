using Gridap

include("../src/initialise.jl")


function area_test(manifold_model,order::Int)
  computed_area = get_surface_area(manifold_model,order)
  _true_area = true_area(get_manifold_name(manifold_model))
  area_diff = abs(computed_area - _true_area)

  println("\nTrue area: $(_true_area)\nComputed area: $(computed_area)\nDifference: $(area_diff)")

  @test area_diff < 1e-8
end

## test the area of cube  -- apply 1 level of refinement
_manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(1),cube)
manifold_model = Adaptivity.refine(_manifold_model)

area_test(manifold_model,2)
area_test(manifold_model,4)


# ## test the area of cubed sphere -- apply 1 level of refinement
# _manifold_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
# manifold_model = Adaptivity.refine(_manifold_model)
# area_test(manifold_model,2)
# area_test(manifold_model,4)
# area_test(manifold_model,6)
# area_test(manifold_model,8)
