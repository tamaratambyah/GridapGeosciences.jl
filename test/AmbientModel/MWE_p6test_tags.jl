using GridapGeosciences
using Gridap
using GridapDistributed
using PartitionedArrays
using MPI
using GridapP4est
using FillArrays
using Gridap.Geometry

dir = @__DIR__
MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
ranks = distribute_with_mpi(LinearIndices((prod(nprocs),)))


n_ref_lvls = 2
radius = 1.0
thickness = 0.19

tags = ["top_boundary"]

#### Parametric model: All good
panel_omodel = CubedSphere3DParametricOctreeDistributedDiscreteModel(
    ranks,radius,thickness;
    num_horizontal_uniform_refinements=n_ref_lvls,
    num_vertical_uniform_refinements=0)
panel_model = panel_omodel.parametric_dmodel
Γ_panel = BoundaryTriangulation(panel_model,tags=tags)

cell_geo_map = geo_map_func(get_forward_map_generator(panel_model),get_panel_ids(Γ_panel))
writevtk_with_cell_geomap(cell_geo_map,Γ_panel,dir*"/boundary",append=false)

### Convert to ambient model
ambient_dmodel = CubedSphereAmbientDistributedDiscreteModel(panel_omodel.parametric_dmodel)

################################################################################
model = ambient_dmodel.models.item

ambient_btrian = BoundaryTriangulation(model,tags=tags)

cmap = get_cell_map(ambient_btrian)
pts = get_cell_ref_coordinates(ambient_btrian)
lazy_map(evaluate,cmap,pts)


writevtk(ambient_btrian,dir*"/ambient_boundary",append=false)
