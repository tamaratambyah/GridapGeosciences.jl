using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays
import Gridap.TensorValues: meas

a = 1.0

include("coarse_cube_surface_2D.jl")
include("Adaptivity/panel_ids_from_refinement.jl")

include("Geometry/Bump.jl")
include("Geometry/PanelRotation.jl")

include("helpers.jl")
include("ManifoldGrid.jl")
include("ManifoldDiscreteModel.jl")

dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


cube_model_2D = UnstructuredDiscreteModel(coarse_cube_surface_2D(a)...)

manifold_grid = ManifoldGrid(cube_model_2D,cube)
panel_ids = get_panel_ids(manifold_grid)

get_cell_node_ids(manifold_grid)
get_cell_permutations(get_grid_topology(cube_model_2D),0)

test_cell_maps(get_cell_map(manifold_grid),get_cell_ref_coordinates(manifold_grid),get_cell_coordinates(manifold_grid))

writevtk(manifold_grid,dir*"/grid",append=false)


manifold_model = ManifoldDiscreteModel(cube_model_2D,cube)

ambient_model = get_ambient_model(manifold_model)
Ω_ambient = Triangulation(ambient_model)
cf_alpha = CellField(x->x[1],Ω_ambient)
cf_beta = CellField(x->x[2],Ω_ambient)
writevtk(Ω_ambient,dir*"/ambient",cellfields=["a"=>cf_alpha,"b"=>cf_beta],append=false)




#### debug
Dc = num_cell_dims(cube_model_2D)
Dp = num_point_dims(cube_model_2D)
panel_ids = get_panel_ids(cube_model_2D)

cube_grid = get_grid(cube_model_2D)
topo = get_grid_topology(cube_model_2D)
cmaps = get_cell_map(cube_grid)
cell_node_ids = get_cell_node_ids(cube_grid)
node_coordinates = get_node_coordinates(cube_grid)
cell_coords = get_cell_coordinates(cube_model_2D)
ref_coords = get_cell_ref_coordinates(cube_grid)

### make the parametric grid
parametric_cell_map = lazy_map(m->cmaps[1],cmaps)
parametric_cell_coords = lazy_map(m->cell_coords[1],panel_ids)

test_cell_maps(parametric_cell_map,ref_coords,parametric_cell_coords)

parametric_grid = UnstructuredGrid(node_coordinates,cell_node_ids,
    get_reffes(cube_grid),get_cell_type(cube_grid),OrientationStyle(cube_grid),
    nothing,parametric_cell_map)

parametric_model = UnstructuredDiscreteModel(parametric_grid,topo,get_face_labeling(cube_model_2D))

num_cells(cube_grid)

## make ambient grid
g =  BumpField(A_bump,B_bump,b_bump)
k = map(p-> PanelRotationField(r1p[p]) ∘ g, panel_ids)

ambient_cell_map = lazy_map(∘,k,parametric_cell_map)

parametric_cell_coords_3D = lazy_map(BumpMap(), parametric_cell_coords)
ambient_cell_coords = lazy_map(R1pPanelMap(), parametric_cell_coords_3D, panel_ids)
ambient_nodes = get_nodes_from_coords(cube_grid,ambient_cell_coords)

test_cell_maps(ambient_cell_map,ref_coords,ambient_cell_coords)


ambient_grid = UnstructuredGrid(ambient_nodes,cell_node_ids,
    get_reffes(cube_grid),get_cell_type(cube_grid),OrientationStyle(cube_grid),
    nothing,ambient_cell_map)

ambient_topo = UnstructuredGridTopology(ambient_grid)
ambient_model = UnstructuredDiscreteModel(ambient_grid,ambient_topo,get_face_labeling(cube_model_2D))

test_cell_maps(ambient_model)





##### test refinement
method = Adaptivity.RedGreenRefinement()
cells_to_refine = nothing

model = cube_model_2D
ctopo = get_grid_topology(model)
coarse_labels = get_face_labeling(model)
# Create new model
rrules, faces_list = Adaptivity.setup_edge_based_rrules(method, model.grid_topology,cells_to_refine)
topo   = Adaptivity.refine_edge_based_topology(ctopo,rrules,faces_list)
reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
topo_grid   = UnstructuredGrid(
  get_vertex_coordinates(topo),
  get_faces(topo,Dc,0),reffes,
  get_cell_type(topo),
  OrientationStyle(topo)
)
glue = Adaptivity.blocked_refinement_glue(rrules)
labels = Adaptivity.refine_face_labeling(coarse_labels,glue,model.grid_topology,topo)

panel_ids = get_panel_ids(model)
new_to_old_cells = glue.n2o_faces_map[Dc+1]
ref_panel_ids = panel_ids[new_to_old_cells]

ref_cell_coords = get_cell_coordinates(topo_grid)
ref_cell_coords[ref_panel_ids.==1]


ref_model = UnstructuredDiscreteModel(topo_grid,topo,labels)
ref_model = AdaptedDiscreteModel(ref_model,model,glue)
_ref_model = AdaptedDiscreteModel(ref_model.model,model,ref_model.glue)

### next level of refinment

ref_model = Adaptivity.refine(cube_model_2D)
ref_cell_coords = get_cell_coordinates(ref_model)
ref_panel_ids = get_panel_ids(ref_model)
ref_cell_coords[ref_panel_ids.==1]


ref_ref_model = Adaptivity.refine(ref_model)
ref_ref_cell_coords = get_cell_coordinates(ref_ref_model)
ref_ref_panel_ids = get_panel_ids(ref_ref_model)
ref_ref_cell_coords[ref_ref_panel_ids.==1]
