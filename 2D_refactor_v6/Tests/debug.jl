αβ = Point(1.0,1.0)
X = forward_map(αβ,1)

_forward_map(p) = αβ -> forward_map(αβ,p)

gradient(_forward_map(1))(αβ)
forward_jacobian(αβ,1)


####################### refine debug

cube_model = coarse_cube_model(π/4,6)
panel_model = parametric_model(cube_model)


#### enter refine routine
model = panel_model
cells_to_refine=nothing
method = Adaptivity.RedGreenRefinement()


## function
panel_ids = get_panel_ids(panel_model)

panel_grid = get_grid(panel_model)
_ctopo = get_grid_topology(panel_model)
Dc = num_cell_dims(_ctopo)

### make the cube topology with proper 3D nodes
panel_cmaps = get_cell_map(panel_grid) ### must get cmaps from the grid NOT THE MODEL!
panel2cube_map = lazy_map(p-> MyAffineField(A_panel2cube[p],b_panel2cube[p]), panel_ids)
cube_cmaps = lazy_map(∘,panel2cube_map,panel_cmaps)

ref_points = get_cell_ref_coordinates(panel_grid)
c_cube_coords = lazy_map(evaluate,cube_cmaps,ref_points)
c_cube_nodes = get_nodes_from_coords(_ctopo,c_cube_coords)

ctopo = UnstructuredGridTopology(c_cube_nodes,
  get_faces(_ctopo,Dc,0),get_cell_type(_ctopo),get_polytopes(_ctopo),OrientationStyle(_ctopo))

c_cube_grid   = Geometry.UnstructuredGrid(
    c_cube_nodes,get_cell_node_ids(panel_grid),get_reffes(panel_grid),
    get_cell_type(panel_grid),OrientationStyle(panel_grid),nothing,cube_cmaps)

writevtk(c_cube_grid,dir*"/cube_grid",append=false)


### apply refinement to the cube
# Create new model
rrules, faces_list = Adaptivity.setup_edge_based_rrules(method, ctopo,cells_to_refine)
topo   = Adaptivity.refine_edge_based_topology(ctopo,rrules,faces_list)
reffes = Adaptivity.map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
ref_cube_grid   = UnstructuredGrid(
  get_vertex_coordinates(topo),
  get_faces(topo,Dc,0),reffes,
  get_cell_type(topo),
  OrientationStyle(topo)
)
glue = Adaptivity.blocked_refinement_glue(rrules)

ref_model = UnstructuredDiscreteModel(ref_cube_grid,topo,FaceLabeling(topo))
acube_model =  AdaptedDiscreteModel(ref_model,cube_model,glue)

writevtk(ref_cube_grid,dir*"/cube_grid",append=false)

writevtk(acube_model,dir*"/cube_model",append=false)


###### make again the panel grid

#### need ref_panel ids
new_to_old_cells = glue.n2o_faces_map[Dc+1]
ref_panel_ids = panel_ids[new_to_old_cells]


#### need to make panel cmaps again
ref_cube_cmaps = get_cell_map(ref_cube_grid)
cube2panel_map = lazy_map(p->MatMultField( A_cube2panel[p] ), ref_panel_ids)
ref_panel_cmaps = lazy_map(∘,cube2panel_map,ref_cube_cmaps)

ref_cube_nodes = get_node_coordinates(ref_cube_grid)

## create the panel model
## to correctly trigger Dc=2,Dp=2, need to have 2D nodes.
## the panel_nodes are just junk 2D nodes, never used by Gridap
## the panel_grid has the bespoke panel_cmap
## the panel_topo is the same as the cube, but with the 2D nodes so that Dp=2
## the panel_labels are from the panel_topo
ref_panel_nodes = map(x->Point(x[2],x[3]),ref_cube_nodes) # these are just junk nodes, never used
ref_panel_grid = Geometry.UnstructuredGrid(ref_panel_nodes,get_cell_node_ids(ref_cube_grid),
          get_reffes(ref_cube_grid),get_cell_type(ref_cube_grid),OrientationStyle(ref_cube_grid),
                    nothing,ref_panel_cmaps)
ref_panel_topo = UnstructuredGridTopology(ref_panel_nodes,get_cell_node_ids(ref_cube_grid),
      get_cell_type(topo),get_polytopes(topo),OrientationStyle(topo))
ref_panel_labels = FaceLabeling(ref_panel_topo)


ref_model = UnstructuredDiscreteModel(ref_panel_grid,ref_panel_topo,ref_panel_labels)
amodel = AdaptedDiscreteModel(ref_model,model,glue)

get_panel_ids(amodel)
