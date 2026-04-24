"""
In this module, we test that the surface area computed by the distributed serial model
is equivalent to the surface area computed by the 2D P4est model at the same
level of refinement on more than 1 processor
i.e. surface area = ∫ᵧ 1 = ∫ 1 √g
"""

module DistributedSurfaceAreaTests

using Gridap
using GridapGeosciences
using GridapP4est
using Test



function compute_surface_area(model, degree::Int)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)

  meas_cf = panelwise_cellfield(sqrtg,Ω)
  surface_area = sum( ∫( 1.0*meas_cf )dΩ )
  return surface_area
end

function main(distribute,nprocs)
  ranks = distribute(LinearIndices((nprocs,)))

  n_ref_lvls = 3
  for radius in [1,2]
    dist_models = get_distributed_refined_models(ranks,nprocs,n_ref_lvls,radius)
    p4est_models = get_octree_refined_models(ranks,n_ref_lvls,radius)
    for degree in collect([2,4,6,8])
      for (d_model,p4_model) in zip(dist_models,p4est_models)
        radius = get_radius(d_model)
        extact_area = 4*π*radius^2

        ### d_model
        d_area = compute_surface_area(d_model, degree)

        ### p4est model
        p4_area = compute_surface_area(p4_model, degree)

        @test d_area ≈ p4_area

        e = abs(d_area-extact_area)/extact_area
        @test e < 1e-2
      end
    end
  end

  @test true
end


end # module
