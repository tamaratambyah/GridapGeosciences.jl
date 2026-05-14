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


n_ref_lvls = 1
radius = 1.0
thickness = 0.19

tags = ["bottom_boundary"]

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

lmodel = ambient_dmodel.models.item
labels = get_face_labeling(lmodel)
labels.tag_to_name
D = num_cell_dims(lmodel)
face_to_mask = Gridap.Geometry.get_face_mask(labels,tags,D-1)

face_to_bgface = findall(face_to_mask)
bgface_to_lcell = Fill(1,num_facets(lmodel))
Γ_ambient = BoundaryTriangulation(lmodel,face_to_bgface,bgface_to_lcell)

topo = get_grid_topology(lmodel)
bgface_grid = Grid(ReferenceFE{D-1},lmodel)

face_grid = view(bgface_grid,face_to_bgface)
cell_grid = get_grid(lmodel)
glue = Gridap.Geometry.FaceToCellGlue(topo,cell_grid,face_grid,face_to_bgface,bgface_to_lcell)
trian = BodyFittedTriangulation(lmodel,face_grid,face_to_bgface)
num_cells(trian)
writevtk(trian,dir*"/boundary_ambient",append=false)
#### nothing there! But there is 24 cells in the trian

cmap = get_cell_map(trian)
pts = get_cell_ref_coordinates(trian)
lazy_map(evaluate,cmap,pts)
#### why is the evalation of the cmap 0??
