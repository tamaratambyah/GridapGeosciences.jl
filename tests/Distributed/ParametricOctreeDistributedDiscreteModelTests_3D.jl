using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed

using DrWatson

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

dir = datadir("Omodel_3D")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

num_horizontal_uniform_refinements = 3
num_vertical_uniform_refinements = 2
o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
                                           num_horizontal_uniform_refinements=num_horizontal_uniform_refinements,
                                           num_vertical_uniform_refinements=num_vertical_uniform_refinements);


panel_model = o3model.parametric_dmodel
Ω_panel =  Triangulation(panel_model)
panel_ids = get_panel_ids(Ω_panel)
panel_ids = get_panel_ids(panel_model)

## Try plotting in distributed
cell_geo_map = geo_map_func(Ω_panel)
writevtk(Ω_panel,dir*"/extruded_model",append=false,geo_map=cell_geo_map)
writevtk(Ω_panel,dir*"/extruded_model",append=false)


################################################################################
#### Plot the boundaries
#### tags = ["intermediate_boundary",  "bottom_boundary",  "top_boundary"]
#### Need to call cell_geo_map map with the panel_ids of Γ
################################################################################
tags = ["bottom_boundary"]
Γ = BoundaryTriangulation(panel_model,tags=tags)
cell_geo_map = geo_map_func(get_panel_ids(Γ))
writevtk(Γ,dir*"/boundary_bottom",append=false,geo_map=cell_geo_map)

tags = ["top_boundary"]
Γ = BoundaryTriangulation(panel_model,tags=tags)
Γ.trians.item_ref[].parent.glue.face_to_bgface
cell_geo_map = geo_map_func(get_panel_ids(Γ))
writevtk(Γ,dir*"/boundary_top",append=false,geo_map=cell_geo_map)

tags = ["intermediate_boundary"]
Γ = BoundaryTriangulation(panel_model,tags=tags)
cell_geo_map = geo_map_func(get_panel_ids(Γ))
writevtk(Γ,dir*"/boundary_intermediate",append=false,geo_map=cell_geo_map)


################################################################################
#### Vertical refinement
################################################################################

# level 0
o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
  num_horizontal_uniform_refinements=0, num_vertical_uniform_refinements=0);
panel_model = o3model.parametric_dmodel
num_cells(panel_model)
Ω_panel = Triangulation(panel_model)
get_panel_ids(panel_model)
get_panel_ids(Ω_panel)
cell_geo_map = geo_map_func(Ω_panel)
writevtk(Ω_panel,dir*"/vmodel0",append=false, geo_map=cell_geo_map)

# level 1
o3model = GridapGeosciences.Distributed.vertically_uniformly_refine(o3model)
panel_model = o3model.parametric_dmodel
num_cells(panel_model)
Ω_panel = Triangulation(panel_model)
get_panel_ids(panel_model)
get_panel_ids(Ω_panel)
cell_geo_map = geo_map_func(Ω_panel)
writevtk(Ω_panel,dir*"/vmodel1",append=false, geo_map=cell_geo_map)


# level 2
o3model = GridapGeosciences.Distributed.vertically_uniformly_refine(o3model)
panel_model = o3model.parametric_dmodel
num_cells(panel_model)
Ω_panel = Triangulation(panel_model)
get_panel_ids(panel_model)
get_panel_ids(Ω_panel)
cell_geo_map = geo_map_func(Ω_panel)
writevtk(Ω_panel,dir*"/vmodel2",append=false, geo_map=cell_geo_map)



################################################################################
#### Horizontal refinement
################################################################################
# level 0
o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
  num_horizontal_uniform_refinements=1, num_vertical_uniform_refinements=1);
panel_model = o3model.parametric_dmodel
num_cells(panel_model)
Ω_panel = Triangulation(panel_model)
get_panel_ids(panel_model)
get_panel_ids(Ω_panel)
cell_geo_map = geo_map_func(Ω_panel)
writevtk(Ω_panel,dir*"/hmodel0",append=false, geo_map=cell_geo_map)

# level 1
o3model = GridapGeosciences.Distributed.horizontally_uniformly_refine(o3model)
panel_model = o3model.parametric_dmodel
num_cells(panel_model)
Ω_panel = Triangulation(panel_model)
get_panel_ids(panel_model)
get_panel_ids(Ω_panel)
cell_geo_map = geo_map_func(Ω_panel)
writevtk(Ω_panel,dir*"/hmodel1",append=false, geo_map=cell_geo_map)


# level 2
o3model = GridapGeosciences.Distributed.horizontally_uniformly_refine(o3model)
panel_model = o3model.parametric_dmodel
num_cells(panel_model)
Ω_panel = Triangulation(panel_model)
get_panel_ids(panel_model)
get_panel_ids(Ω_panel)
cell_geo_map = geo_map_func(Ω_panel)
writevtk(Ω_panel,dir*"/hmodel2",append=false, geo_map=cell_geo_map)
