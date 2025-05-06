using Gridap
include("initalise.jl")


manifold_model = ManifoldDiscreteModel(cube_model_2D,cube)
ref_manifold_model = Adaptivity.refine(manifold_model)

@test is_child(ref_manifold_model,manifold_model)

writevtk(get_ambient_grid(get_grid(ref_manifold_model)),dir*"/ref_grid",append=false)

### debugg
model = manifold_model
method = Adaptivity.RedGreenRefinement()
cells_to_refine = nothing

Dc = num_cell_dims(model)
ctopo = get_grid_topology(model)
coarse_labels = get_face_labeling(model)
# Create new model
rrules, faces_list = Adaptivity.setup_edge_based_rrules(method, model.grid_topology,cells_to_refine)
topo   = Adaptivity.refine_edge_based_topology(ctopo,rrules,faces_list)
reffes = map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
cube_grid   = UnstructuredGrid(
  get_vertex_coordinates(topo),
  get_faces(topo,Dc,0),reffes,
  get_cell_type(topo),
  OrientationStyle(topo)
)
glue = Adaptivity.blocked_refinement_glue(rrules)
labels = Adaptivity.refine_face_labeling(coarse_labels,glue,model.grid_topology,topo)

# make the manifold grid
panel_ids = get_panel_ids(model)
new_to_old_cells = glue.n2o_faces_map[Dc+1]
ref_panel_ids = panel_ids[new_to_old_cells]


node_coordinates = get_node_coordinates(cube_grid)
cell_coords = get_cell_coordinates(cube_grid)
cell_node_ids = get_cell_node_ids(cube_grid)
cell_coords[ref_panel_ids.==1]
cell_node_ids[findall(ref_panel_ids.==1)]



panel1_nodes = get_panel_1_nodes_from_coords(cube_grid,cell_coords,ref_panel_ids)
