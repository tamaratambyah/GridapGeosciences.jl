"""
In this module, we test that the surface area computed by the serial model
is equivalent to the surface area computed by the 2D P4est model at the same
level of refinement (on 1 processor)
i.e. surface area = ∫ᵧ 1 = ∫ 1 √g
"""

using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed

include("../convergence_tools.jl")

function compute_surface_area(model, degree::Int)
  Ω = Triangulation(model)
  panel_ids = get_panel_ids(model)
  dΩ = Measure(Ω,degree)

  meas_cf = panelwise_cellfield(sqrtg,Ω,panel_ids)
  surface_area = sum( ∫( 1.0*meas_cf )dΩ )
  return surface_area
end

function main(distribute,nprocs)
  ranks = distribute(LinearIndices((nprocs,)))

  n_ref_lvls = 3
  serial_models = get_refined_models(n_ref_lvls)
  dist_models = get_octree_refined_models(ranks,n_ref_lvls)
  for degree in collect([2,4,6,8])
    for (s_model,d_model) in zip(serial_models,dist_models)
      extact_area = 4*π*RADIUS^2

      ### s_model
      s_area = compute_surface_area(s_model, degree)

      ### d_model
      d_area = compute_surface_area(d_model, degree)

      @test s_area ≈ d_area

      e = abs(s_area-extact_area)/extact_area
      @test e < 1e-2

    end

  end


end

with_mpi() do distribute
  main(distribute,1)
end
