module DistributedPanelIdsTests
using Gridap
using GridapGeosciences
using Gridap.Adaptivity
using Test
using GridapDistributed
# using Test; @testset "DistributedPanelIds" begin include("tests/Distributed/DistributedPanelIdsTests.jl") end

include("../convergence_tools.jl")

################################################################################
## Test the panel ids from the BodyFittedTriangulation trian are the same as the
## panel ids from the model
################################################################################

nprocs = 6
ranks = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end

n_ref_lvls = 2
dmodels = get_distributed_refined_models(ranks,nprocs,n_ref_lvls)

function test_distributed_panel_ids(dpanel_model)
  trian = Triangulation(dpanel_model)
  gids = get_cell_gids(dpanel_model)
  map(get_panel_ids(dpanel_model), get_panel_ids(trian), local_views(dpanel_model)) do p1, p2, model
    @test p1 == p2
    @test length(p1) == num_cells(model)
  end

  map(get_owned_panel_ids(dpanel_model),get_panel_ids(trian),partition(gids)) do p1, p2, cid
    owned_cells = own_to_local(cid)
    _p2 = p2[owned_cells]
    @test p1 == _p2
  end

end


test_distributed_panel_ids(dmodels[1])

end
