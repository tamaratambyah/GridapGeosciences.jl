using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed
using DrWatson

## should return ghost+owned panel ids
MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

dir = datadir("Omodel_2D_refinement")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

# level 0
omodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=0)
panel_model = omodel.parametric_dmodel
num_cells(panel_model)
Ω_panel = Triangulation(panel_model)
get_panel_ids(panel_model)
get_panel_ids(Ω_panel)
cell_geo_map = geo_map_func(Ω_panel)
writevtk(Ω_panel,dir*"/model0",append=false, geo_map=cell_geo_map)

function adapt_model(ranks,model::ParametricOctreeDistributedDiscreteModel)
  cell_partition=get_cell_gids(model.octree_dmodel)
  ref_flags=map(ranks,partition(cell_partition)) do rank,indices
      flags=zeros(Cint,length(indices))
      flags.=refine_flag
  end
  Gridap.Adaptivity.adapt(omodel,ref_flags)
end

# level 1
omodel = adapt_model(ranks,omodel)
panel_model = omodel.parametric_dmodel
num_cells(panel_model)
Ω_panel = Triangulation(panel_model)
get_panel_ids(panel_model)
get_panel_ids(Ω_panel)
cell_geo_map = geo_map_func(Ω_panel)
writevtk(Ω_panel,dir*"/model1",append=false, geo_map=cell_geo_map)


# level 2
omodel = adapt_model(ranks,omodel)
panel_model = omodel.parametric_dmodel
num_cells(panel_model)
Ω_panel = Triangulation(panel_model)
get_panel_ids(panel_model)
get_panel_ids(Ω_panel)
cell_geo_map = geo_map_func( Ω_panel)
writevtk(Ω_panel,dir*"/model2",append=false, geo_map=cell_geo_map)

# level 3
omodel = adapt_model(ranks,omodel)
panel_model = omodel.parametric_dmodel
num_cells(panel_model)
Ω_panel = Triangulation(panel_model)
get_panel_ids(panel_model)
get_panel_ids(Ω_panel)
cell_geo_map = geo_map_func( Ω_panel)
writevtk(Ω_panel,dir*"/model3",append=false, geo_map=cell_geo_map)
