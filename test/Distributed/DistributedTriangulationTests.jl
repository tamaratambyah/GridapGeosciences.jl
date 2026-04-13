"""
Test the construction of meshes by evaluating the cellmaps
"""

module DistributedTriangulationTests

using Gridap
using GridapGeosciences
using GridapDistributed
using GridapP4est
using Test

################################################################################
## Test the evaluation of cmaps on DistributedTriangulations
## i.e. are all the cellmaps there
################################################################################

function test_triangulation(trian::GridapDistributed.DistributedTriangulation)
  map(trian.trians) do trian
    cmap = get_cell_map(trian)
    pts = get_cell_ref_coordinates(trian)
    lazy_map(evaluate,cmap,pts)
    @test true
  end
end

function main(distribute,nprocs)
  test_distributedParametricDiscreteModel(distribute,nprocs)
  test_ParametricOctreeDistributedDiscreteModel(distribute,nprocs)
  test_Parametric3DOctreeDistributedDiscreteModel(distribute,nprocs)
end

function test_distributedParametricDiscreteModel(distribute,nprocs)
  ranks = distribute(LinearIndices((nprocs,)))

  # i_am_main(ranks) && println("--test DistributedParametricDiscreteModel")

  n_ref_lvls = 2
  dmodels = get_distributed_refined_models(ranks,nprocs,n_ref_lvls)

  model = dmodels[2]

  trian = Triangulation(model)
  test_triangulation(trian)

  btrian = BoundaryTriangulation(model)
  test_triangulation(btrian)

  strian = SkeletonTriangulation(model)
  test_triangulation(strian)

  @test true
end


function test_ParametricOctreeDistributedDiscreteModel(distribute,nprocs)
  ranks = distribute(LinearIndices((nprocs,)))

  # i_am_main(ranks) && println("--test ParametricOctreeDistributedDiscreteModel")

  n_ref_lvls = 2
  omodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=n_ref_lvls)
  model = omodel.parametric_dmodel

  trian = Triangulation(model)
  test_triangulation(trian)

  btrian = BoundaryTriangulation(model)
  test_triangulation(btrian)

  strian = SkeletonTriangulation(model)
  test_triangulation(strian)

  @test true
end


function test_Parametric3DOctreeDistributedDiscreteModel(distribute,nprocs)

  ranks = distribute(LinearIndices((nprocs,)))

  # i_am_main(ranks) && println("--test 3D Parametric3DOctreeDistributedDiscreteModel")
  n_ref_lvls = 2
  o3model = Parametric3DOctreeDistributedDiscreteModel(ranks;
  num_horizontal_uniform_refinements=n_ref_lvls, num_vertical_uniform_refinements=n_ref_lvls);
  panel_model = o3model.parametric_dmodel

  trian = Triangulation(panel_model)
  test_triangulation(trian)

  strian = SkeletonTriangulation(panel_model)
  test_triangulation(strian)



  tags = ["bottom_boundary"]
  Γ = BoundaryTriangulation(panel_model,tags=tags)
  cell_geo_map = geo_map_func(get_panel_ids(Γ))
  test_triangulation(Γ)
  # writevtk(Γ,dir*"/boundary_bottom",append=false,geo_map=cell_geo_map)

  tags = ["top_boundary"]
  Γ = BoundaryTriangulation(panel_model,tags=tags)
  Γ.trians.item_ref[].parent.glue.face_to_bgface
  cell_geo_map = geo_map_func(get_panel_ids(Γ))
  test_triangulation(Γ)
  # writevtk(Γ,dir*"/boundary_top",append=false,geo_map=cell_geo_map)

  tags = ["intermediate_boundary"]
  Γ = BoundaryTriangulation(panel_model,tags=tags)
  cell_geo_map = geo_map_func(get_panel_ids(Γ))
  test_triangulation(Γ)
  # writevtk(Γ,dir*"/boundary_intermediate",append=false,geo_map=cell_geo_map)

  @test true
end


end ## module
