using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed

using DrWatson
dir = datadir("Distributed")
!isdir(dir) && mkdir(dir)

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

num_horizontal_uniform_refinements = 1
num_vertical_uniform_refinements = 1
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


##### testing boundary tags
labels = panel_model.models.item.face_labeling
labs = Gridap.Geometry.get_tag_name(labels)



# tags = ["intermediate_boundary",  "bottom_boundary",  "top_boundary"]
tags = ["bottom_boundary"]
Γ = BoundaryTriangulation(panel_model,tags=tags)
Γ.trians.item_ref[].parent.glue.face_to_bgface
cell_geo_map = geo_map_func(Γ)
writevtk(Γ,dir*"/boundary_bottom",append=false,geo_map=cell_geo_map)

tags = ["top_boundary"]
Γ = BoundaryTriangulation(panel_model,tags=tags)
Γ.trians.item_ref[].parent.glue.face_to_bgface
cell_geo_map = geo_map_func(Γ)
writevtk(Γ,dir*"/boundary_top",append=false,geo_map=cell_geo_map)

tags = ["intermediate_boundary"]
Γ = BoundaryTriangulation(panel_model,tags=tags)
cell_geo_map = geo_map_func(Γ)
writevtk(Γ,dir*"/boundary_intermediate",append=false)#,geo_map=cell_geo_map)

o3model = GridapGeosciences.Distributed.vertically_uniformly_refine(o3model)
panel_model = o3model.parametric_dmodel

Ω_panel =  Triangulation(panel_model)
panel_ids = get_panel_ids(Ω_panel)
panel_ids = get_panel_ids(panel_model)

## Try plotting in distributed
cell_geo_map = geo_map_func(Ω_panel)

map(local_views(panel_model)) do model
	println(typeof(model))
    Ω_panel=Triangulation(model)
	cell_geo_map = geo_map_func(Ω_panel)
	writevtk(Ω_panel,dir*"/vrefined_extruded_model",append=false,geo_map=cell_geo_map)
end
# writevtk(Ω_panel,dir*"/vrefined_extruded_model",append=false,geo_map=cell_geo_map)
# writevtk(Ω_panel,dir*"/verfined_extruded_model",append=false)