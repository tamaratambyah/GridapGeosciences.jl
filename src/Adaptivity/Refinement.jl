"""
Refine CubedSphereParametricDiscreteModel

Jordi's refinement routine requires nodes that are properly placed in space.
This is because he determines the location of new nodes as the middle of two
existing nodes.
Thus, must construct the topology to have the proper 3D nodes on the cube.
From there, reconstruct the panel model.

The method is as follows: given a panel model
1. reconstrcut the cube using panel2cube mappings, and extract the 3D nodes
2. reconstrcut the topology using 3D nodes
3. apply EdgeBasedRefinement as per Jordi's routine to obtain the refined_cube_grid
   and adaptivity glue
4. construct refined_panel_grid using refined_cube_grid, and reconstrcut CubedSphereParametricDiscreteModel

"""

import Gridap.Adaptivity: EdgeBasedRefinement

function Gridap.Adaptivity.refine(model::CubedSphereParametricDiscreteModel,args...;refinement_method="red_green",kwargs...)
  return Gridap.Adaptivity.refine(Gridap.Adaptivity.string_to_refinement(refinement_method, model),model,args...;kwargs...)
end

function Gridap.Adaptivity.refine(method::EdgeBasedRefinement,model::CubedSphereParametricDiscreteModel;cells_to_refine=nothing)
  Dc = num_cell_dims(model)
  panel_grid = get_grid(model)
  _ctopo = get_grid_topology(model)
  c_panel_ids = get_panel_ids(model)

  ### make the topology with proper 3D nodes. The 3D nodes are the nodes on the cube
  ### so apply panel2cube mapping
  panel_cmaps = get_cell_map(panel_grid) ### must get cmaps from the grid NOT THE MODEL!
  panel2cube_map = lazy_map(p-> MyAffineField(A_panel2cube[p],b_panel2cube[p]), c_panel_ids)
  cube_cmaps = lazy_map(∘,panel2cube_map,panel_cmaps)

  cell_ref_coords = get_cell_ref_coordinates(panel_grid)
  c_cube_coords = lazy_map(evaluate,cube_cmaps,cell_ref_coords)
  c_cube_nodes = get_nodes_from_coords(_ctopo,c_cube_coords)

  ctopo = UnstructuredGridTopology(c_cube_nodes,
    get_faces(_ctopo,Dc,0),get_cell_type(_ctopo),get_polytopes(_ctopo),OrientationStyle(_ctopo))

  ### Now enter the refinement routine provided by Jordi to obtain:
  ### 1. the refined cubed grid (i.e. refined cmaps)
  ### 2. the glue (to get refined panel ids)
  rrules, faces_list = Gridap.Adaptivity.setup_edge_based_rrules(method, ctopo,cells_to_refine)
  topo   = Gridap.Adaptivity.refine_edge_based_topology(ctopo,rrules,faces_list)
  reffes = Gridap.Adaptivity.map(p->LagrangianRefFE(Float64,p,1),get_polytopes(topo))
  r_cube_grid   = UnstructuredGrid(
    get_vertex_coordinates(topo),
    get_faces(topo,Dc,0),reffes,
    get_cell_type(topo),
    OrientationStyle(topo)
  )
  glue = Gridap.Adaptivity.blocked_refinement_glue(rrules)

  ###### Construct the refined panel grid
  ### get the r_panel ids, and make r_panel_cmaps
  new_to_old_cells = glue.n2o_faces_map[Dc+1]
  r_panel_ids = c_panel_ids[new_to_old_cells]

  r_cube_cmaps = get_cell_map(r_cube_grid)
  cube2panel_map = lazy_map(p->MatMultField( A_cube2panel[p] ), r_panel_ids)
  r_panel_cmaps = lazy_map(∘,cube2panel_map,r_cube_cmaps)

  ## create the refined panel model
  ## to correctly trigger Dc=2,Dp=2, need to have 2D nodes.
  ## the panel_nodes are just junk 2D nodes, never used by Gridap
  ## the panel_grid has the bespoke panel_cmap
  ## the panel_topo is the same as the cube, but with the 2D nodes so that Dp=2
  ## the panel_labels are from the panel_topo
  r_cube_nodes = get_node_coordinates(r_cube_grid)
  r_panel_nodes = map(x->Point(x[2],x[3]),r_cube_nodes) # these are just junk nodes, never used
  r_panel_grid = Gridap.Geometry.UnstructuredGrid(r_panel_nodes,get_cell_node_ids(r_cube_grid),
            get_reffes(r_cube_grid),get_cell_type(r_cube_grid),OrientationStyle(r_cube_grid),
                      nothing,r_panel_cmaps)
  r_panel_topo = UnstructuredGridTopology(r_panel_nodes,get_cell_node_ids(r_cube_grid),
        get_cell_type(topo),get_polytopes(topo),OrientationStyle(topo))
  r_panel_labels = FaceLabeling(r_panel_topo)


  ### construct the refined model and return adapted model
  ref_model = CubedSphere2DParametricDiscreteModel(r_panel_grid,r_panel_topo,r_panel_labels,r_panel_ids,get_radius(model))
  return AdaptedDiscreteModel(ref_model,model,glue)
end


### return true if the number of cells is the same, since === will not return
### true in the case of CubedSphereParametricDiscreteModel
function Gridap.Adaptivity.is_child(m1::AdaptedDiscreteModel,m2::CubedSphereParametricDiscreteModel)
  # println("parametric is child")
  return num_cells(get_parent(m1)) === num_cells(m2)
end



"""
Refine CubedSphereAmbientDiscreteModel

Refine the associated CubedSphereParametricDiscreteModel as above
"""

function Gridap.Adaptivity.refine(model::CubedSphereAmbientDiscreteModel,args...;refinement_method="red_green",kwargs...)
  panel_model = get_parametric_model(model)
  ref_panel_model = Gridap.Adaptivity.refine(panel_model)
  CubedSphereAmbientDiscreteModel(ref_panel_model)
end
