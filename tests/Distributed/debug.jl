using Gridap
using GridapGeosciences
using DrWatson
using FillArrays
using Gridap.Geometry, Gridap.CellData, Gridap.Fields, Gridap.Helpers
using Test


using GridapDistributed
using PartitionedArrays
using MPIPreferences
MPIPreferences.use_jll_binary()


dir = datadir("Distributed")
include("../convergence_tools.jl")


n_ref_lvls = 2

nprocs = 6

ranks = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end

dmodels = get_distributed_refined_models(ranks,nprocs,n_ref_lvls)



################################################################################
#### BodyFittedTriangulation
################################################################################


dpanel_model = dmodels[2]
get_panel_ids(dpanel_model)

trian = Triangulation(dpanel_model)

map(btrian.trians) do trian
  cmap = get_cell_map(trian)
  pts = get_cell_ref_coordinates(trian)
  lazy_map(evaluate,cmap,pts)
  @test true
end

cell_geo_map = geo_map_func(get_owned_panel_ids(trian))
writevtk(trian,dir*"/distributed_model",append=false,geo_map=cell_geo_map)

################################################################################
#### BoundaryTriangulation
#### Need to return the mask that is all the interior cells
################################################################################

btrian = BoundaryTriangulation(dpanel_model)
b_panel_ids =  get_panel_ids(btrian)

b_geo_map = geo_map_func(b_panel_ids)
writevtk(btrian,dir*"/distributed_boundary_trian",append=false,geo_map=b_geo_map)



################################################################################
#### SkeletonTriangulation
#### Need to dispatch to serial
################################################################################
skel = SkeletonTriangulation(dpanel_model)

skel_panel_ids = get_panel_ids(skel)
_skel_panel_ids = get_skel_panel_ids(skel_panel_ids)

skel_geo_map = geo_map_func(_skel_panel_ids)
writevtk(skel,dir*"/distributed_skel_trian",append=false,geo_map=skel_geo_map)
