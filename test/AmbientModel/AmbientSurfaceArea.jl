"""
In this module, we test that the surface area computed by the
   serial parametric model: = ∫ᵧ 1 = ∫ 1 √g
is equvialent to the surface area computed by the
   serial ambient model: = ∫ᵧ 1
"""

module AmbientSurfaceArea

using Gridap
using GridapGeosciences
using GridapP4est
using Test

import GridapGeosciences.Geometry: ParametricModels

function compute_surface_area(
  ambient_model::Union{CubedSphereAmbientDiscreteModel,CubedSphereAmbientDistributedDiscreteModel{2,3,<:CubedSphereAmbientDiscreteModel}},
  degree::Int)
  Ω = Triangulation(ambient_model)
  dΩ = Measure(Ω,degree)

  surface_area = sum( ∫( 1.0 )dΩ )
  return surface_area
end

function compute_surface_area(
  model::Union{ParametricModels,CubedSphere2DParametricDistributedDiscreteModel},
  degree::Int)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  meas_cf = ParametricCellField(sqrtg,Ω)
  surface_area = sum( ∫( 1.0*meas_cf )dΩ )
  return surface_area
end


function main(ambient_models::AbstractArray;_i_am_main=true)
  for degree in collect([2,4,6,8])
    for (ambient_model) in ambient_models
      panel_model = get_parametric_model(ambient_model)
      radius = get_radius(panel_model)
      extact_area = 4*π*radius^2

      ### panel_model
      panel_area = compute_surface_area(panel_model, degree)
      e_panel = abs(panel_area-extact_area)/extact_area
      _i_am_main && println("Parametric error:", e_panel)
      @test e_panel < 1e-2

      ### ambient_model
      ambient_area = compute_surface_area(ambient_model, degree)
      e_ambient = abs(ambient_area-extact_area)/extact_area
      _i_am_main && println("Ambient error:", e_ambient)
      @test e_ambient < 1e-2

      e_comparison = e_ambient - e_panel
      _i_am_main && println("Comparison error:", e_comparison)
      @test e_comparison < 1e-12
    end
  end

end


end # module
