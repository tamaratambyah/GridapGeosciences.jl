"""
In this module, we test that the surface area computed by the
   serial parametric model: = ∫ᵧ 1 = ∫ 1 √g
is equvialent to the surface area computed by the
   serial ambient model: = ∫ᵧ 1
"""

module SurfaceArea

using Gridap
using GridapGeosciences
using GridapP4est
using Test

import GridapGeosciences.Geometry: ParametricModels

function compute_surface_area(model::CubedSphereAmbientDiscreteModel, degree::Int)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  surface_area = sum( ∫( 1.0 )dΩ )
  return surface_area
end

function compute_surface_area(model::ParametricModels, degree::Int)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  meas_cf = ParametricCellField(sqrtg,Ω)
  surface_area = sum( ∫( 1.0*meas_cf )dΩ )
  return surface_area
end


function main(serial_panel_models::AbstractArray{<:ParametricModels})
  for degree in collect([2,4,6,8])
    for (s_panel_model) in serial_panel_models
      radius = get_radius(s_panel_model)
      extact_area = 4*π*radius^2

      ### s_panel_model
      s_panel_area = compute_surface_area(s_panel_model, degree)
      e_panel = abs(s_panel_area-extact_area)/extact_area
      println("Parametric error:", e_panel)
      @test e_panel < 1e-2

      ### s_ambient_model
      s_ambient_model = CubedSphereAmbientDiscreteModel(s_panel_model)
      s_ambient_area = compute_surface_area(s_ambient_model, degree)
      e_ambient = abs(s_ambient_area-extact_area)/extact_area
      println("Ambient error:", e_ambient)
      @test e_ambient < 1e-2

      e_comparison = e_ambient - e_panel
      println("Comparison error:", e_comparison)
      @test e_comparison < 1e-12
    end
  end

end



end # module
