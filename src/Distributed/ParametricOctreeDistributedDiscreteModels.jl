"""
	ParametricOctreeDistributedDiscreteModel{Dc,Dp}
    
	generic_dmodel includes a distributed array of proc-local ParametricDiscreteModel{Dc,Dp} objects

	octree_model manages the topology of the parametric mesh 
"""

const NPANELS = 6


# Ideas:
# * GenericDistributedDiscreteModel should be created out of OctreeDistributedDiscreteModel
# * Many methods of ParametricOctreeDistributedDiscreteModel should be forwarded to the underlying
#   GenericDistributedDiscreteModel 
# * At some point, we should create the 3D coordinates of the nodes in the cube surface 
# * To this end, we can use a cell-wise, DG-like array, to be transformed into a nodal-like array,
#   generated out of the composition of two maps:
#    - map from the reference cell to the local coordinate system of the panel [-1,1]^2 or [0,1]^2
#    - map from the local coordinate system of the panel to the 3D cube surface 

struct ParametricOctreeDistributedDiscreteModel{2,2,
	                                            A<:OctreeDistributedDiscreteModel{2,2},
												B<:GenericDistributedDiscreteModel{2,2}} <: DistributedDiscreteModel{2,2}
  a::Float64
  octree_dmodel::A
  parametric_dmodel::B
end 

# TO-THINK: not sure if radius of the cubed sphere is required?
# Right now, in the serial version of the code, we are passing the length of the cube edges
# (see coarse_cube_surface_3D function)               
function ParametricOctreeDistributedDiscreteModel(ranks; a=π/4, num_initial_uniform_refinements=0)
	coarse_model = _create_parametric_octree_dmodel_coarse_model()
	octree_dmodel, cell_wise_vertex_2D_coordinates, cell_panels =
	        _generate_octree_dmodel_2D_coordinates_and_panels(ranks, coarse_model, num_initial_uniform_refinements);
    
	# Transform 2D coordinates to 3D coordinates on the cube surface
	# ...

	# Build the proc-local ParametricDiscreteModels
	# ...

	# Build the GenericDistributedDiscreteModel
	# ...

	# ParametricOctreeDistributedDiscreteModel(a, octree_dmodel)
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

  function setup_2D_to_3D_coase_cell_map(a)
	 nodes_3d = _CCAM_cube_nodes_3d(a)
	 cell_node_ids = _CCAM_panel_wise_node_ids(NPANELS)
	 cell_vertices_3D = setup_coarse_cell_vertices_3D_coordinates(cell_node_ids,nodes_3d)
     scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(QUAD,Gridap.ReferenceFEs.lagrangian,Float64,1)
	 cell_shape_funs = FillArrays.Fill( Gridap.Fields.get_shape_functions(scalar_reffe), NPANELS) 
	 lazy_map(linear_combination,cell_to_coords,cell_to_shapefuns)
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


# Dc will be always 2, anyway
function _generate_octree_dmodel_2D_coordinates_and_panels(ranks, 
	                                                       coarse_model::DiscreteModel{Dc}, 
														   num_initial_uniform_refinements,
														   coarse_cell_wise_vertex_coordinates,
														   coarse_cell_panel) where Dc
   comm = ranks.comm
   pXest_type = GridapP4est._dim_to_pXest_type(Dc)

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
   cell_coordinates_and_panels=generate_cell_coordinates_and_panels(ranks,
                                             coarse_discrete_model,
                                             coarse_cell_wise_vertex_coordinates,
                                             coarse_cell_panel,
                                             ptr_pXest_connectivity,
                                             ptr_pXest,
                                             ptr_pXest_ghost)	
   cube_grid_top=generate_cube_grid_top(cell_corner_lids_nlcorners)
   models= map(cube_grid_top) do cube_grid_top
            Gridap.Geometry.UnstructuredDiscreteModel(cube_grid_top)
         end
   dmodel=GridapDistributed.DistributedDiscreteModel(models,cellindices)
   dmodel, cell_coordinates_and_panels...
end

