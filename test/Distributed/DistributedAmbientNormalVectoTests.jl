module DistributedAmbientNormalTests

using Gridap
using GridapDistributed
using GridapGeosciences
using Test


function test_debug_equiv(p,q)
  map(p,q  )do p,q
    dif = p - q
    max_dif = map(x->maximum(norm.(x)),dif)
    @test all(max_dif .< 1e-12)
  end
end


function main(distribute,nprocs)

  ranks = distribute(LinearIndices((nprocs,)))

  n_ref_lvls = 2
  radius = 1.0
  ambient_model = CubedSphereAmbientDistributedDiscreteModel(
    ranks, radius; num_initial_uniform_refinements=n_ref_lvls)
  panel_model = get_parametric_model(ambient_model)


  ##############################################################################
  ########## Ambient model
  ##############################################################################
  Ω_ambient = Triangulation(ambient_model)
  n_surface_ambient = get_surface_normal(Ω_ambient)

  Λ_ambient = SkeletonTriangulation(ambient_model)
  n_Λ_ambient = get_normal_vector(Λ_ambient)
  pts_ambient = get_cell_points(Λ_ambient)

  ### Test the skeleton normal is in the tangent space of the sphere
  normal_component = (n_Λ_ambient⋅n_surface_ambient).plus(pts_ambient)
  max_dif = map(x->maximum(norm.(x)),normal_component)
  map(max_dif) do _d
    @test all(_d .< 1e-12)
  end


  normal_component = (n_Λ_ambient⋅n_surface_ambient).minus(pts_ambient)
  max_dif = map(x->maximum(norm.(x)),normal_component)
  map(max_dif) do _d
    @test all(_d .< 1e-12)
  end

  ##############################################################################
  ########## Parametric model
  ##############################################################################

  Λ_panel = SkeletonTriangulation(panel_model)
  panel_ids = get_panel_ids(panel_model)
  forward_map_generator = get_forward_map_generator(panel_model)
  cell_geo_map = geo_map_func(forward_map_generator,panel_ids)
  n_Λ_mapped = pushforward_normal(Λ_panel,cell_geo_map)
  pts_panel = get_cell_points(Λ_panel)

  ## Test equivalence with ambient model
  p = n_Λ_mapped.plus(pts_panel)
  q = n_Λ_ambient.plus(pts_ambient)
  test_debug_equiv(p,q)

  p = n_Λ_mapped.minus(pts_panel)
  q = n_Λ_ambient.minus(pts_ambient)
  test_debug_equiv(p,q)



  @test true

end

end # module
