using Gridap
include("src/initialise.jl")




a = π/4
r = 1.0*sqrt(3.0)
global A_bump, B_bump, b_bump = bump_matrics(a)

cube_grid_3D,topo_3D,face_labels_3D,panel_ids = coarse_cube_surface_3D(a)


Dc = num_cell_dims(cube_grid_3D)
Dp_amb = num_point_dims(cube_grid_3D)

cube_cell_coords_3D = get_cell_coordinates(cube_grid_3D)

coords_panel1_3D = lazy_map(Rp1PanelMap3D(), cube_cell_coords_3D, panel_ids)
coords_panel1_2D = lazy_map(BumpMap(), coords_panel1_3D)

latlon_panel1 = lazy_map(GnomonicMap(), coords_panel1_2D)
sphere_panel1 = lazy_map(SigmaMap(r),latlon_panel1)
sphere_panelp = lazy_map(R1pPanelMap3D(), sphere_panel1, panel_ids)



ambient_cell_coords = lazy_map(RoundVectorValues(),sphere_panelp)

_cube_grid_3D, = coarse_cube_surface_3D(1.0)
get_cell_coordinates(_cube_grid_3D)


ambient_nodes = get_nodes_from_coords(cube_grid_3D,ambient_cell_coords)

ambient_grid = Gridap.Geometry.UnstructuredGrid(ambient_nodes,get_cell_node_ids(cube_grid_3D),
    get_reffes(cube_grid_3D),get_cell_type(cube_grid_3D),OrientationStyle(cube_grid_3D))

cmaps = get_cell_map(ambient_grid)
evaluate(cmaps[1],Point(1,1))



manifold_grid = ManifoldGrid(cube,cube_grid_3D,panel_ids)

num_point_dims(manifold_grid)
num_cell_dims(manifold_grid)


model = ManifoldDiscreteModel(manifold_grid)
num_point_dims(model)

ref_model = Adaptivity.refine(model)
num_point_dims(ref_model)
num_cell_dims(ref_model)

writevtk(get_ambient_grid(get_grid(ref_model)),dir*"/ref_grid",append=false)

get_panel_ids(ref_model)

ref_ref_model = Adaptivity.refine(ref_model)
num_point_dims(ref_ref_model)
num_cell_dims(ref_ref_model)
get_panel_ids(ref_ref_model)
writevtk(get_ambient_grid(get_grid(ref_ref_model)),dir*"/ref_ref_grid",append=false)
