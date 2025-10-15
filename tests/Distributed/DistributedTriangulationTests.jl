"""
Test the construction of meshes by evaluating the cellmaps
"""

module DistributedTriangulationTests
using Gridap
using GridapGeosciences
using Test
using GridapDistributed

include("../convergence_tools.jl")

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
  ranks = distribute(LinearIndices((nprocs,)))

  n_ref_lvls = 2
  dmodels = get_distributed_refined_models(ranks,nprocs,n_ref_lvls)

  model = dmodels[2]

  trian = Triangulation(model)
  test_triangulation(trian)

  btrian = BoundaryTriangulation(model)
  test_triangulation(btrian)

  strian = SkeletonTriangulation(model)
  test_triangulation(strian)
end


end ## module
