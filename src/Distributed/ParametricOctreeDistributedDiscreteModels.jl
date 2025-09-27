"""
	ParametricOctreeDistributedDiscreteModel{Dc,Dp}
    
	generic_dmodel includes a distributed array of proc-local ParametricDiscreteModel{Dc,Dp} objects

	octree_model manages the topology of the parametric mesh 
"""

const NPANELS = 6
const CUBE_SURFACE_HALF_EDGE = π/4 


# Ideas:
# * GenericDistributedDiscreteModel should be created out of OctreeDistributedDiscreteModel
# * Many methods of ParametricOctreeDistributedDiscreteModel should be forwarded to the underlying
#   GenericDistributedDiscreteModel 
# * At some point, we should create the 3D coordinates of the nodes in the cube surface 
# * To this end, we can use a cell-wise, DG-like array, to be transformed into a nodal-like array,
#   generated out of the composition of two maps:
#    - map from the reference cell to the local coordinate system of the panel [-1,1]^2 or [0,1]^2
#    - map from the local coordinate system of the panel to the 3D cube surface 

struct ParametricOctreeDistributedDiscreteModel{A<:OctreeDistributedDiscreteModel{2,2}, 
	                                            B<:GenericDistributedDiscreteModel{2,2}} <: DistributedDiscreteModel{2,2}   
												
  octree_dmodel::A
  parametric_dmodel::B
end 

# TO-THINK: not sure if radius of the cubed sphere is required?
# Right now, in the serial version of the code, we are passing the length of the cube edges
# (see coarse_cube_surface_3D function)               
function ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=0)
	coarse_model = _create_parametric_octree_dmodel_coarse_model()
	octree_dmodel, cell_wise_vertex_2D_coordinates, cell_panels =
	        _generate_octree_dmodel_2D_coordinates_and_panels(ranks, 
			                                                  coarse_model, 
															  num_initial_uniform_refinements, 
															  setup_coarse_cell_vertices_2D_coordinates(), 
															  collect(1:NPANELS))

	# Transform 2D coordinates to 3D coordinates on the cube surface
	map_2D_to_3D_panels = setup_2D_to_3D_coarse_cell_map()
    map_2D_to_3D_cells  = map(cell_panels) do cell_panels
	   lazy_map(Reindex(map_2D_to_3D_panels), cell_panels)
	end
    cell_wise_vertex_3D_coordinates = 
	  map(cell_wise_vertex_2D_coordinates, map_2D_to_3D_cells) do cell_wise_vertex_2D_coordinates,
			                                                      map_2D_to_3D_cells
	   lazy_map(evaluate, map_2D_to_3D_cells, cell_wise_vertex_2D_coordinates)
	  end


	# Build the proc-local ParametricDiscreteModels
	parametric_models = map(local_views(octree_dmodel.dmodel), 
	                        cell_wise_vertex_3D_coordinates,
							cell_panels) do omodel, cell_wise_vertex_3D_coordinates, cell_panels

        cube_cmaps = setup_2D_to_3D_cell_map(cell_wise_vertex_3D_coordinates)
		panel_cmaps = setup_panel_cmaps(cube_cmaps, cell_panels)

        ogrid = get_grid(omodel)
		otopo = get_grid_topology(omodel)
		panel_grid = UnstructuredGrid(get_node_coordinates(ogrid),
		                                       get_cell_node_ids(ogrid),
											   get_reffes(ogrid),
											   get_cell_type(ogrid),
											   OrientationStyle(ogrid),
                                               nothing,
											   panel_cmaps)
        panel_topo = UnstructuredGridTopology(get_node_coordinates(ogrid),
		                                      get_cell_node_ids(ogrid),
											  get_cell_type(ogrid),
											  get_polytopes(otopo),
											  OrientationStyle(ogrid))
        panel_labels = FaceLabeling(panel_topo)

        ParametricDiscreteModel(panel_grid,
		                        panel_topo,
								panel_labels,
								cell_panels)                                                  
	end

	# Build the GenericDistributedDiscreteModel
	generic_dmodel = GenericDistributedDiscreteModel(parametric_models, get_cell_gids(octree_dmodel.dmodel))
	ParametricOctreeDistributedDiscreteModel(octree_dmodel, generic_dmodel)
end 

function _create_parametric_octree_dmodel_coarse_model()
  # 6 panels (cells), 4 corners (vertices) each panel
  cell_node_ids = _CCAM_panel_wise_node_ids(NPANELS)
  node_coordinates = Vector{Point{2,Float64}}(undef,8)
  # These coordinates are junk coordinates, but they have to be unique.
  # This a precondition for the OctreeDistributedDiscreteModel to properly
  # count the vertices of the coarse model
  for i in 1:length(node_coordinates)
    node_coordinates[i]=Point{2,Float64}(Float64(i),Float64(i))
  end
  polytope=QUAD
  scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
  cell_types=collect(Fill(1,length(cell_node_ids)))
  cell_reffes=[scalar_reffe]
  grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                          cell_node_ids,
                                          cell_reffes,
                                          cell_types,
                                          Gridap.Geometry.NonOriented())
  Gridap.Geometry.UnstructuredDiscreteModel(grid)
end



  # This info goes to P4est ...
  function setup_coarse_cell_vertices_2D_coordinates()
	# We use the reference coordinate system of the QUAD
	# to build the 2D->3D map for each panel of the cube surface 
    data = vcat([ QUAD.vertex_coords for i=1:6 ]...)
    ptr = [1, 5, 9, 13, 17, 21, 25]
    Gridap.Arrays.Table(data,ptr)
  end 

  # This info is required for the parametric model's map
  function setup_coarse_cell_vertices_3D_coordinates(cell_node_ids, nodes_3d)
     nodes_3d_cell_wise = [ nodes_3d[cell_node_ids.data[i]]  for i=1:length(cell_node_ids.data) ]
     ptr = [1, 5, 9, 13, 17, 21, 25]
     Gridap.Arrays.Table(nodes_3d_cell_wise,ptr)
  end

function setup_2D_to_3D_coarse_cell_map()
  nodes_3d = _CCAM_cube_nodes_3d(CUBE_SURFACE_HALF_EDGE)
  cell_node_ids = _CCAM_panel_wise_node_ids(NPANELS)
  cell_vertices_3D = setup_coarse_cell_vertices_3D_coordinates(cell_node_ids,nodes_3d)
  scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(QUAD,Gridap.ReferenceFEs.lagrangian,Float64,1)
  cell_shape_funs = FillArrays.Fill( Gridap.ReferenceFEs.get_shapefuns(scalar_reffe), NPANELS) 
  lazy_map(linear_combination,cell_vertices_3D,cell_shape_funs)
end

function setup_2D_to_3D_cell_map(cell_vertices_3D)
  scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(QUAD,Gridap.ReferenceFEs.lagrangian,Float64,1)
  cell_shape_funs = FillArrays.Fill( Gridap.ReferenceFEs.get_shapefuns(scalar_reffe), length(cell_vertices_3D)) 
  lazy_map(linear_combination,cell_vertices_3D,cell_shape_funs)
end

function generate_cube_grid_top(cell_vertex_lids_nlvertices)
  map(cell_vertex_lids_nlvertices[1],cell_vertex_lids_nlvertices[2]) do cell_vertex_lids,nlvector
     node_coordinates=Vector{Point{2,Float64}}(undef,nlvector)
     polytope=QUAD
     scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
     cell_types=collect(Fill(1,length(cell_vertex_lids)))
     cell_reffes=[scalar_reffe]
     cell_vertex_lids_gridap=Gridap.Arrays.Table(cell_vertex_lids.data,cell_vertex_lids.ptrs)
     grid=Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                      cell_vertex_lids_gridap,
                                      cell_reffes,
                                      cell_types,
                                      Gridap.Geometry.NonOriented())
     grid
  end
end

function set_coarse_cell_vertices_coordinates!( pconn :: Ptr{P4est_wrapper.p4est_connectivity_t},
                                                coarse_discrete_model :: DiscreteModel{2,2},
                                                panel,
                                                ref_cell_coordinates)
  @assert panel ≤ num_cells(coarse_discrete_model)
  @assert panel ≥ 1
  trian=Triangulation(coarse_discrete_model)
  cell_vertices=Gridap.Geometry.get_cell_node_ids(trian)
  #println(cell_vertices)
  cell_vertices_panel=cell_vertices[panel]
  conn=pconn[]
  vertices=unsafe_wrap(Array,
                       conn.vertices,
                       length(Gridap.Geometry.get_node_coordinates(coarse_discrete_model))*3)
  for (l,g) in enumerate(cell_vertices_panel)
     vertices[(g-1)*3+1]=ref_cell_coordinates[l][1]
     vertices[(g-1)*3+2]=ref_cell_coordinates[l][2]
  end
end

function generate_cell_coordinates_and_panels(parts,
                                   coarse_discrete_model,
                                   coarse_cell_wise_vertex_coordinates,
                                   coarse_cell_panel,
                                   ptr_pXest_connectivity,
                                   ptr_pXest,
                                   ptr_pXest_ghost)

  Dc=2
  PXEST_CORNERS=4
  pXest_ghost = ptr_pXest_ghost[]
  pXest = ptr_pXest[]

  # Obtain ghost quadrants
  ptr_ghost_quadrants = Ptr{P4est_wrapper.p4est_quadrant_t}(pXest_ghost.ghosts.array)

  tree_offsets = unsafe_wrap(Array, pXest_ghost.tree_offsets, pXest_ghost.num_trees+1)
  dcell_coordinates_and_panels=map(parts) do part
     ncells=pXest.local_num_quadrants+pXest_ghost.ghosts.elem_count
     panels = Vector{Int}(undef,ncells)
     data = Vector{Point{Dc,Float64}}(undef,ncells*PXEST_CORNERS)
     ptr  = generate_ptr(ncells)
     current=1
     current_cell=1
     vxy=Vector{Cdouble}(undef,Dc)
     pvxy=pointer(vxy,1)
     for itree=1:pXest_ghost.num_trees
       tree = p4est_tree_array_index(pXest.trees, itree-1)[]
       if tree.quadrants.elem_count > 0
          set_coarse_cell_vertices_coordinates!( ptr_pXest_connectivity, 
                                                 coarse_discrete_model, 
                                                 itree, 
                                                 coarse_cell_wise_vertex_coordinates[itree])
       end
       for cell=1:tree.quadrants.elem_count
          panels[current_cell]=coarse_cell_panel[itree]
          quadrant=p4est_quadrant_array_index(tree.quadrants, cell-1)[]
          for vertex=1:PXEST_CORNERS
            GridapP4est.p4est_get_quadrant_vertex_coordinates(ptr_pXest_connectivity,
                                                  p4est_topidx_t(itree-1),
                                                  quadrant.x,
                                                  quadrant.y,
                                                  quadrant.level,
                                                  Cint(vertex-1),
                                                  pvxy)
            data[current]=Point{Dc,Float64}(vxy...)
            current=current+1
          end
          current_cell=current_cell+1
       end
     end

     # Go over ghost cells
     for i=1:pXest_ghost.num_trees
      if tree_offsets[i+1]-tree_offsets[i] > 0
        set_coarse_cell_vertices_coordinates!( ptr_pXest_connectivity, 
                                               coarse_discrete_model, 
                                               i, 
                                               coarse_cell_wise_vertex_coordinates[i])
      end
      for j=tree_offsets[i]:tree_offsets[i+1]-1
          panels[current_cell]=i
          quadrant = ptr_ghost_quadrants[j+1]
          for vertex=1:PXEST_CORNERS
            GridapP4est.p4est_get_quadrant_vertex_coordinates(ptr_pXest_connectivity,
                                                     p4est_topidx_t(i-1),
                                                     quadrant.x,
                                                     quadrant.y,
                                                     quadrant.level,
                                                     Cint(vertex-1),
                                                     pvxy)

          #  if (MPI.Comm_rank(comm.comm)==0)
          #     println(vxy)
          #  end
          data[current]=Point{Dc,Float64}(vxy...)
          current=current+1
         end
         current_cell=current_cell+1
       end
     end
     Gridap.Arrays.Table(data,ptr), panels
  end |> tuple_of_arrays
end


function _generate_octree_dmodel_2D_coordinates_and_panels(ranks, 
	                                                       coarse_model::DiscreteModel{2,2}, 
														   num_uniform_refinements,
														   coarse_cell_wise_vertex_coordinates,
														   coarse_cell_panel)
   comm = ranks.comm
   pXest_type = GridapP4est._dim_to_pXest_type(2)

   ptr_pXest_connectivity,
      ptr_pXest,
        ptr_pXest_ghost,
          ptr_pXest_lnodes = GridapP4est.setup_ptr_pXest_objects(pXest_type,
                                                                 comm,
                                                                 coarse_model,
                                                                 num_uniform_refinements)

   cellindices = GridapP4est.setup_cell_prange(pXest_type,ranks,ptr_pXest,ptr_pXest_ghost)
   cell_vertex_gids=GridapP4est.generate_cell_vertex_gids(ptr_pXest_lnodes,cellindices)
   cell_corner_lids=GridapP4est.generate_cell_corner_lids(cell_vertex_gids)
   cell_corner_lids_nlcorners=map(cell_corner_lids) do cell_corner_lids
     cell_corner_lids,maximum(cell_corner_lids.data)
   end |> tuple_of_arrays
   cell_coordinates, panels=generate_cell_coordinates_and_panels(ranks,
                                             coarse_model,
                                             coarse_cell_wise_vertex_coordinates,
                                             coarse_cell_panel,
                                             ptr_pXest_connectivity,
                                             ptr_pXest,
                                             ptr_pXest_ghost)

	GridapP4est.pXest_lnodes_destroy(pXest_type,ptr_pXest_lnodes)
    GridapP4est.pXest_ghost_destroy(pXest_type,ptr_pXest_ghost)
	
   cube_grid_top=generate_cube_grid_top(cell_corner_lids_nlcorners)
   models= map(cube_grid_top) do cube_grid_top
            Gridap.Geometry.UnstructuredDiscreteModel(cube_grid_top)
         end
   dmodel=GridapDistributed.DistributedDiscreteModel(models,cellindices)

   non_conforming_glue = GridapP4est._create_conforming_model_non_conforming_glue(dmodel)


   omodel=OctreeDistributedDiscreteModel(2,
                                  2,
                                  ranks,
                                  dmodel,
                                  non_conforming_glue,
                                  coarse_model,
                                  ptr_pXest_connectivity,
                                  ptr_pXest,
                                  pXest_type,
                                  GridapP4est.PXestUniformRefinementRuleType(),
                                  true,
                                  nothing)


   omodel, cell_coordinates, panels
end

