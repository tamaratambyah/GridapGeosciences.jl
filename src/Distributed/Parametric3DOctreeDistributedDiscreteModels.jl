struct Parametric3DOctreeDistributedDiscreteModel{A<:OctreeDistributedDiscreteModel{3,3}, 
                                                  B<:GenericDistributedDiscreteModel{3,3}} <: DistributedDiscreteModel{3,3}   
                                                
  octree_dmodel::A
  parametric_dmodel::B
end 
              
function Parametric3DOctreeDistributedDiscreteModel(ranks; 
                                                    num_horizontal_uniform_refinements=0,
                                                    num_vertical_uniform_refinements=0)


    
    coarse_model = _create_parametric_octree_dmodel_coarse_model()

    octree_dmodel, cell_wise_vertex_alpha_beta_gamma_coordinates, cell_panels =
            _generate_octree_dmodel_alpha_beta_gamma_coordinates_and_panels(ranks, 
                                                                coarse_model, 
                                                                num_horizontal_uniform_refinements, 
                                                                num_vertical_uniform_refinements,
                                                                setup_coarse_cell_vertices_alpha_beta_coordinates(), 
                                                                collect(1:NPANELS))

    # Build the proc-local ParametricDiscreteModels
    parametric_models = _setup_parametric_models(octree_dmodel, 
                                                 cell_wise_vertex_alpha_beta_gamma_coordinates,
                                                  cell_panels)

    # Build the GenericDistributedDiscreteModel
    generic_dmodel = GenericDistributedDiscreteModel(parametric_models, get_cell_gids(octree_dmodel.dmodel))
    Parametric3DOctreeDistributedDiscreteModel(octree_dmodel, generic_dmodel)
end 

function _setup_parametric_models(octree_dmodel::OctreeDistributedDiscreteModel{3,3}, 
                                 cell_wise_vertex_alpha_beta_gamma_coordinates,
                                 cell_panels)

    map(local_views(octree_dmodel.dmodel), 
                            cell_wise_vertex_alpha_beta_gamma_coordinates,
                            cell_panels) do omodel, cell_wise_vertex_alpha_beta_gamma_coordinates, cell_panels

        alpha_beta_cmap = setup_alpha_beta_gamma_cell_map(cell_wise_vertex_alpha_beta_gamma_coordinates)

        ogrid = get_grid(omodel)
        otopo = get_grid_topology(omodel)
        panel_grid = UnstructuredGrid(get_node_coordinates(ogrid),
                                               get_cell_node_ids(ogrid),
                                               get_reffes(ogrid),
                                               get_cell_type(ogrid),
                                               OrientationStyle(ogrid),
                                               nothing,
                                               alpha_beta_cmap)
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
end 



function setup_alpha_beta_gamma_cell_map(cell_vertices_alpha_beta_gamma)
  scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(HEX,Gridap.ReferenceFEs.lagrangian,Float64,1)
  cell_shape_funs = 
     FillArrays.Fill( Gridap.ReferenceFEs.get_shapefuns(scalar_reffe), length(cell_vertices_alpha_beta_gamma) )
  lazy_map(linear_combination,cell_vertices_alpha_beta_gamma,cell_shape_funs)
end

function generate_3D_cube_grid_top(cell_vertex_lids_nlvertices)
  map(cell_vertex_lids_nlvertices[1],cell_vertex_lids_nlvertices[2]) do cell_vertex_lids,nlvector
     node_coordinates=Vector{Point{3,Float64}}(undef,nlvector)
     polytope=HEX
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

function generate_cell_alpha_beta_gamma_coordinates_and_panels(parts,
                                   coarse_discrete_model,
                                   coarse_coarse_cell_wise_vertex_alpha_beta_coordinates,
                                   coarse_cell_panel,
                                   ptr_pXest_connectivity,
                                   ptr_pXest,
                                   ptr_pXest_ghost)

  Dc=3
  PXEST_CORNERS=2^Dc
  pXest_ghost = ptr_pXest_ghost[]
  pXest = ptr_pXest[]
  pXest_type = GridapP4est.P6estType()

  dcell_coordinates_and_panels=map(parts) do part
     panels = Int[]
     data = Point{Dc,Float64}[]
     ncells = 0
     vxy=Vector{Cdouble}(undef,Dc)
     pvxy=pointer(vxy,1)
     for itree=1:pXest_ghost.num_trees
       tree = GridapP4est.pXest_tree_array_index(pXest_type, pXest, itree-1)[]
       if tree.quadrants.elem_count > 0
          set_coarse_cell_vertices_coordinates!( ptr_pXest_connectivity[].conn4, 
                                                 coarse_discrete_model, 
                                                 itree, 
                                                 coarse_coarse_cell_wise_vertex_alpha_beta_coordinates[itree])
       end
       for cell=1:tree.quadrants.elem_count
          quadrant=GridapP4est.pXest_quadrant_array_index(pXest_type,tree,cell-1)[]
          
          # Loop over layers in the current column
          for l=1:GridapP4est.pXest_num_quadrant_layers(pXest_type,quadrant)
            layer=GridapP4est.pXest_get_layer(pXest_type, quadrant, pXest, l-1)
            coords=GridapP4est.pXest_cell_coords(pXest_type,quadrant,layer)
            levels=GridapP4est.pXest_get_quadrant_and_layer_levels(pXest_type,quadrant,layer)
            push!(panels, coarse_cell_panel[itree])
            for vertex=1:PXEST_CORNERS
              GridapP4est.pXest_get_quadrant_vertex_coordinates(pXest_type,
                                                               ptr_pXest_connectivity,
                                                               p4est_topidx_t(itree-1),
                                                               coords,
                                                               levels,
                                                               Cint(vertex-1),
                                                               pvxy)
              push!(data, Point{Dc,Float64}(vxy[3],vxy[1],vxy[2]))
            end
            ncells = ncells+1
          end 
       end
     end

     function sc_array_p4est_locidx_t_index(sc_array_object::sc_array_t, it)
      @assert sc_array_object.elem_size == sizeof(p4est_locidx_t)
      @assert it in 0:sc_array_object.elem_count
      ptr=Ptr{p4est_locidx_t}(sc_array_object.array + sc_array_object.elem_size*it)
      return unsafe_wrap(Array, ptr, 1)[]
     end

     column_ghost = pXest_ghost.column_ghost[]
     ptr_p2est_ghost_quadrants = GridapP4est._unwrap_ghost_quadrants(pXest_type, pXest_ghost)
     ptr_p4est_ghost_quadrants = GridapP4est._unwrap_ghost_quadrants(GridapP4est.P4estType(), column_ghost)

     tree_offsets = unsafe_wrap(Array, column_ghost.tree_offsets, pXest_ghost.num_trees+1)
     
     current_ghost_column=0

     # Go over ghost cells
     for i=1:pXest_ghost.num_trees
       for j=tree_offsets[i]:tree_offsets[i+1]-1
          p4est_quadrant = ptr_p4est_ghost_quadrants[j+1]
          k = sc_array_p4est_locidx_t_index(pXest_ghost.column_layer_offsets[],current_ghost_column)
          l = sc_array_p4est_locidx_t_index(pXest_ghost.column_layer_offsets[],current_ghost_column+1)
          for m=k:l-1
            p2est_quadrant = ptr_p2est_ghost_quadrants[m+1]
            coords=pXest_cell_coords(pXest_type,p4est_quadrant,p2est_quadrant)
            levels=pXest_get_quadrant_and_layer_levels(pXest_type,p4est_quadrant,p2est_quadrant)
            for vertex=1:PXEST_CORNERS
                GridapP4est.pXest_get_quadrant_vertex_coordinates(pXest_type,
                                                      ptr_pXest_connectivity,
                                                      p4est_topidx_t(i-1),
                                                      coords,
                                                      levels,
                                                      Cint(vertex-1),
                                                      pvxy)
                push!(data, Point{Dc,Float64}(vxy[3],vxy[1],vxy[2]))
            end
            ncells=ncells+1
          end
          current_ghost_column=current_ghost_column+1
       end
    end
    ptr=generate_ptr(Dc,ncells)
    Gridap.Arrays.Table(data,ptr), panels
  end |> tuple_of_arrays
end


function _generate_octree_dmodel_alpha_beta_gamma_coordinates_and_panels(ranks, 
                                                                         coarse_model::DiscreteModel{2,2}, 
                                                                         num_horizontal_uniform_refinements,
                                                                         num_vertical_uniform_refinements,
                                                                         coarse_cell_wise_vertex_alpha_beta_coordinates,
                                                                         coarse_cell_panel)
   
   comm = ranks.comm
   Dc=3
   pXest_type = GridapP4est.P6estType()
   pXest_refinement_rule_type = GridapP4est.PXestHorizontalRefinementRuleType()

   extrusion_vector::Vector{Float64}=[0.0,0.0,1.0]
   
   ptr_pXest_connectivity=GridapP4est.setup_pXest_connectivity(pXest_type,
                                                   coarse_model,
                                                   extrusion_vector)
   
   ptr_pXest=GridapP4est.setup_pXest(pXest_type, 
                                    comm, 
                                    ptr_pXest_connectivity, 
                                    num_horizontal_uniform_refinements,
                                    num_vertical_uniform_refinements)


    ptr_pXest_ghost=GridapP4est.setup_pXest_ghost(pXest_type,ptr_pXest)
    ptr_pXest_lnodes=GridapP4est.setup_pXest_lnodes_nonconforming(pXest_type, ptr_pXest, ptr_pXest_ghost)

    cell_prange = GridapP4est.setup_cell_prange(pXest_type, ranks, ptr_pXest, ptr_pXest_ghost)

    gridap_cell_faces,
        non_conforming_glue=
        GridapP4est.generate_cell_faces_and_non_conforming_glue(pXest_type, 
                                                                pXest_refinement_rule_type,
                                                                ptr_pXest_lnodes, 
                                                                cell_prange)
    
    GridapP4est.pXest_lnodes_destroy(pXest_type,ptr_pXest_lnodes)

     nlvertices = map(non_conforming_glue) do ncglue
          ncglue.num_regular_faces[1]+ncglue.num_hanging_faces[1]
     end

     cell_coordinates, panels=generate_cell_alpha_beta_gamma_coordinates_and_panels(ranks,
                                              coarse_model,
                                              setup_coarse_cell_vertices_alpha_beta_coordinates(),
                                              coarse_cell_panel,
                                              ptr_pXest_connectivity,
                                              ptr_pXest,
                                              ptr_pXest_ghost)

    cell_corner_coordinates =
       _generate_zero_cell_corner_coordinates(pXest_type,gridap_cell_faces[1])


    function JaggedToTable(x::MPIArray{<:JaggedArray})
      map(x) do x
        Gridap.Arrays.Table(x.data,x.ptrs)
      end
    end

    grid,topology=_generate_topology_grid_and_topology(pXest_type,
                                                       JaggedToTable(gridap_cell_faces[1]),
                                                       cell_corner_coordinates)

    map(topology,gridap_cell_faces[Dc]) do topology,cell_faces
      cell_faces_gridap = Gridap.Arrays.Table(cell_faces.data,cell_faces.ptrs)
      topology.n_m_to_nface_to_mfaces[Dc+1,Dc] = cell_faces_gridap
      topology.n_m_to_nface_to_mfaces[Dc,Dc+1] = Gridap.Geometry.generate_cells_around(cell_faces_gridap)
    end

    
    map(topology,gridap_cell_faces[Dc-1]) do topology,cell_edges
      cell_edges_gridap = Gridap.Arrays.Table(cell_edges.data,cell_edges.ptrs)
      topology.n_m_to_nface_to_mfaces[Dc+1,Dc-1] = cell_edges_gridap
      topology.n_m_to_nface_to_mfaces[Dc-1,Dc+1] = Gridap.Geometry.generate_cells_around(cell_edges_gridap)
    end 

     models=map(grid, topology) do grid, topology
        labeling=FaceLabeling(topology)
        Gridap.Geometry.UnstructuredDiscreteModel(grid, topology, labeling)
     end
     dmodel=GridapDistributed.DistributedDiscreteModel(models,cell_prange)


     GridapP4est.pXest_ghost_destroy(pXest_type,ptr_pXest_ghost)

    omodel= GridapP4est.OctreeDistributedDiscreteModel(Dc,
                                                       Dc,
                                                       ranks,
                                                       dmodel,
                                                       non_conforming_glue,
                                                       coarse_model,
                                                       ptr_pXest_connectivity,
                                                       ptr_pXest,
                                                       pXest_type,
                                                       pXest_refinement_rule_type,
                                                       true,
                                                       nothing)

    omodel, cell_coordinates, panels
end


