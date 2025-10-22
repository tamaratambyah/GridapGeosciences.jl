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

include("forward_map_3D.jl")


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))
coarse_model = GridapGeosciences.Distributed._create_parametric_octree_dmodel_coarse_model()
add_tag!(coarse_model.face_labeling, "interior", [1])
num_horizontal_uniform_refinements = 0
num_vertical_uniform_refinements = 0

dmodel=AnisotropicallyAdapted3DDistributedDiscreteModel(ranks,
                                                        coarse_model,
														num_horizontal_uniform_refinements,
														num_vertical_uniform_refinements)

dmodel.dmodel.models.item_ref[].grid.node_coordinates
dmodel.dmodel.models.item_ref[].grid_topology.n_m_to_nface_to_mfaces[4,1]

model = dmodel.dmodel
Ω_panel =  Triangulation(model)

panel_ids = get_panel_ids(Ω_panel)
# panel_ids = get_panel_ids(model) ### BROKEN - Not parametric model

### plot parametric model
map(local_views(model)) do model
  writevtk(Triangulation(model),dir*"/model_serial",append=false)
end

## Try plotting in distributed
cell_geo_map = geo_map_func_3D(Ω_panel)
writevtk(Ω_panel,dir*"/extruded_model",append=false,geo_map=cell_geo_map)


## Try plotting in serial
map(local_views(model),panel_ids) do model, panel_ids
    cell_geo_map = lazy_map(p ->  ForwardMapPanel3D(p), panel_ids)
	writevtk(Triangulation(model),dir*"/extruded_model_serial",append=false,geo_map=cell_geo_map)
end
