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
            _generate_octree_alpha_beta_gamma_coordinates_and_panels(ranks,
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

        alpha_beta_gamma_cmap = setup_alpha_beta_gamma_cell_map(cell_wise_vertex_alpha_beta_gamma_coordinates)

        ogrid = get_grid(omodel)
        otopo = get_grid_topology(omodel)
        olabels = Gridap.Geometry.get_face_labeling(omodel)
        panel_grid = UnstructuredGrid(get_node_coordinates(ogrid),
                                     get_cell_node_ids(ogrid),
                                     get_reffes(ogrid),
                                     get_cell_type(ogrid),
                                     OrientationStyle(ogrid),
                                     nothing,
                                     alpha_beta_gamma_cmap)
        ParametricDiscreteModel(panel_grid,
                                otopo,
                                olabels,
                                cell_panels)
    end
end



function setup_alpha_beta_gamma_cell_map(cell_vertices_alpha_beta_gamma)
  scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(HEX,Gridap.ReferenceFEs.lagrangian,Float64,1)
  cell_shape_funs =
     FillArrays.Fill( Gridap.ReferenceFEs.get_shapefuns(scalar_reffe), length(cell_vertices_alpha_beta_gamma) )
  lazy_map(linear_combination,cell_vertices_alpha_beta_gamma,cell_shape_funs)
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
       if tree_offsets[i+1]-tree_offsets[i] > 0
          set_coarse_cell_vertices_coordinates!( ptr_pXest_connectivity[].conn4,
                                                 coarse_discrete_model,
                                                 i,
                                                 coarse_coarse_cell_wise_vertex_alpha_beta_coordinates[i])
       end

       for j=tree_offsets[i]:tree_offsets[i+1]-1
          p4est_quadrant = ptr_p4est_ghost_quadrants[j+1]
          k = sc_array_p4est_locidx_t_index(pXest_ghost.column_layer_offsets[],current_ghost_column)
          l = sc_array_p4est_locidx_t_index(pXest_ghost.column_layer_offsets[],current_ghost_column+1)
          for m=k:l-1
            p2est_quadrant = ptr_p2est_ghost_quadrants[m+1]
            coords=GridapP4est.pXest_cell_coords(pXest_type,p4est_quadrant,p2est_quadrant)
            levels=GridapP4est.pXest_get_quadrant_and_layer_levels(pXest_type,p4est_quadrant,p2est_quadrant)
            push!(panels, coarse_cell_panel[i])
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

function dummy_grid_and_topology_function(pXest_type::GridapP4est.P6estType,
                                          non_conforming_glue,
                                          cell_vertices,
                                          ptr_pXest_connectivity,
                                          ptr_pXest,
                                          ptr_pXest_ghost)
  function JaggedToTable(x::MPIArray{<:JaggedArray})
      map(x) do x
        Gridap.Arrays.Table(x.data,x.ptrs)
      end
  end
  grid,topology=_generate_topology_grid_and_topology(pXest_type,JaggedToTable(cell_vertices))
end


function _generate_octree_alpha_beta_gamma_coordinates_and_panels(ranks,
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

   ptr_pXest=P4est_wrapper.p6est_new_ext(comm,
                  ptr_pXest_connectivity,
                  Cint(0),
                  Cint(num_horizontal_uniform_refinements), # min_level
                  Cint(num_vertical_uniform_refinements),   # min_zlevel
                  Cint(1),                       # num_zroot
                  Cint(1),                       # fill_uniform
                  Cint(1),                       # data_size
                  C_NULL,                        # init_fn
                  C_NULL)                        # user_pointer


    ptr_pXest_ghost=GridapP4est.setup_pXest_ghost(pXest_type,ptr_pXest)
    ptr_pXest_lnodes=GridapP4est.setup_pXest_lnodes_nonconforming(pXest_type, ptr_pXest, ptr_pXest_ghost)


    dmodel,non_conforming_glue  = GridapP4est.setup_non_conforming_distributed_discrete_model(pXest_type,
                                                    GridapP4est.PXestHorizontalRefinementRuleType(),
                                                    ranks,
                                                    coarse_model,
                                                    ptr_pXest_connectivity,
                                                    ptr_pXest,
                                                    ptr_pXest_ghost,
                                                    ptr_pXest_lnodes;
                                                    grid_and_topology_function=dummy_grid_and_topology_function,
                                                    grid_and_topology_bottom_function=dummy_grid_and_topology_function)

    cell_coordinates, panels=generate_cell_alpha_beta_gamma_coordinates_and_panels(ranks,
                                          coarse_model,
                                          setup_coarse_cell_vertices_alpha_beta_coordinates(),
                                          coarse_cell_panel,
                                          ptr_pXest_connectivity,
                                          ptr_pXest,
                                          ptr_pXest_ghost)

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

    GridapP4est.pXest_lnodes_destroy(pXest_type,ptr_pXest_lnodes)
    GridapP4est.pXest_ghost_destroy(pXest_type,ptr_pXest_ghost)

    omodel, cell_coordinates, panels
end


function vertically_uniformly_refine(parametric_model::Parametric3DOctreeDistributedDiscreteModel)
  ptr_new_pXest = GridapP4est._vertically_uniformly_refine!(parametric_model.octree_dmodel)

  pXest_type = parametric_model.octree_dmodel.pXest_type

  # Extract ghost and lnodes
  ptr_pXest_ghost  = GridapP4est.setup_pXest_ghost(pXest_type, ptr_new_pXest)
  ptr_pXest_lnodes = GridapP4est.setup_pXest_lnodes_nonconforming(pXest_type, ptr_new_pXest, ptr_pXest_ghost)

  pXest_refinement_rule_type = parametric_model.octree_dmodel.pXest_refinement_rule_type
  ranks = parametric_model.octree_dmodel.parts
  coarse_model = parametric_model.octree_dmodel.coarse_model
  ptr_pXest_connectivity = parametric_model.octree_dmodel.ptr_pXest_connectivity

  fmodel,non_conforming_glue  = GridapP4est.setup_non_conforming_distributed_discrete_model(pXest_type,
                                              parametric_model.octree_dmodel.pXest_refinement_rule_type,
                                              ranks,
                                              coarse_model,
                                              ptr_pXest_connectivity,
                                              ptr_new_pXest,
                                              ptr_pXest_ghost,
                                              ptr_pXest_lnodes;
                                              grid_and_topology_function=dummy_grid_and_topology_function,
                                              grid_and_topology_bottom_function=dummy_grid_and_topology_function)

  cell_coordinates, panels=generate_cell_alpha_beta_gamma_coordinates_and_panels(ranks,
                                          coarse_model,
                                          setup_coarse_cell_vertices_alpha_beta_coordinates(),
                                          collect(1:NPANELS),
                                          ptr_pXest_connectivity,
                                          ptr_new_pXest,
                                          ptr_pXest_ghost)

   GridapP4est.pXest_ghost_destroy(pXest_type,ptr_pXest_ghost)
   GridapP4est.pXest_lnodes_destroy(pXest_type,ptr_pXest_lnodes)

   pXest_refinement_rule_type = GridapP4est.PXestVerticalRefinementRuleType()
   _refinement_and_coarsening_flags = map(partition(get_cell_gids(parametric_model.octree_dmodel))) do indices
     flags  = Vector{Cint}(undef,length(local_to_global(indices)))
     flags .= refine_flag
   end

   stride = GridapP4est.pXest_stride_among_children(pXest_type,
                                                    pXest_refinement_rule_type,
                                                    parametric_model.octree_dmodel.ptr_pXest)
   adaptivity_glue = GridapP4est._compute_fine_to_coarse_model_glue(pXest_type,
                                                       pXest_refinement_rule_type,
                                                       ranks,
                                                       parametric_model.octree_dmodel.dmodel,
                                                       fmodel,
                                                       _refinement_and_coarsening_flags,
                                                       stride)
   adaptive_models = map(local_views(parametric_model.octree_dmodel),
                           local_views(fmodel),
                           adaptivity_glue) do model, fmodel, glue
           Gridap.Adaptivity.AdaptedDiscreteModel(fmodel,model,glue)
   end
   fmodel = GridapDistributed.GenericDistributedDiscreteModel(adaptive_models,get_cell_gids(fmodel))
   ref_model = OctreeDistributedDiscreteModel(3,3,
                                              ranks,
                                              fmodel,
                                              non_conforming_glue,
                                              coarse_model,
                                              ptr_pXest_connectivity,
                                              ptr_new_pXest,
                                              pXest_type,
                                              parametric_model.octree_dmodel.pXest_refinement_rule_type,
                                              false,
                                              parametric_model)

   # Build the proc-local ParametricDiscreteModels
   parametric_models = _setup_parametric_models(ref_model,
                                                cell_coordinates,
                                                panels)

   adaptive_models = map(parametric_models,
                         local_views(parametric_model.parametric_dmodel),
                         local_views(ref_model.dmodel)) do parametric_model,
                                                           parametric_model_parent,
                                                           octree_dmodel_adapted_model
        Gridap.Adaptivity.AdaptedDiscreteModel(parametric_model,
                                               parametric_model_parent,
                                               get_adaptivity_glue(octree_dmodel_adapted_model))
   end
   generic_dmodel =
      GenericDistributedDiscreteModel(adaptive_models, get_cell_gids(ref_model.dmodel))

   Parametric3DOctreeDistributedDiscreteModel(ref_model, generic_dmodel)
end

function horizontally_uniformly_refine(parametric_model::Parametric3DOctreeDistributedDiscreteModel)

    num_cols = GridapP4est.num_locally_owned_columns(parametric_model.octree_dmodel)
    _refinement_and_coarsening_flags = map(num_cols) do num_cols
        flags  = Vector{Cint}(undef,num_cols)
        flags .= refine_flag
    end
    ptr_new_pXest = GridapP4est._horizontally_refine_coarsen_balance!(parametric_model.octree_dmodel,
                                                                      _refinement_and_coarsening_flags)


    pXest_type = parametric_model.octree_dmodel.pXest_type
    pXest_refinement_rule_type = parametric_model.octree_dmodel.pXest_refinement_rule_type
    ranks = parametric_model.octree_dmodel.parts
    coarse_model = parametric_model.octree_dmodel.coarse_model
    ptr_pXest_connectivity = parametric_model.octree_dmodel.ptr_pXest_connectivity

    # Extract ghost and lnodes
    ptr_pXest_ghost  = GridapP4est.setup_pXest_ghost(pXest_type, ptr_new_pXest)
    ptr_pXest_lnodes = GridapP4est.setup_pXest_lnodes_nonconforming(pXest_type, ptr_new_pXest, ptr_pXest_ghost)

    fmodel,non_conforming_glue  = GridapP4est.setup_non_conforming_distributed_discrete_model(pXest_type,
                                              parametric_model.octree_dmodel.pXest_refinement_rule_type,
                                              ranks,
                                              coarse_model,
                                              ptr_pXest_connectivity,
                                              ptr_new_pXest,
                                              ptr_pXest_ghost,
                                              ptr_pXest_lnodes;
                                              grid_and_topology_function=dummy_grid_and_topology_function,
                                              grid_and_topology_bottom_function=dummy_grid_and_topology_function)

    cell_coordinates, panels=generate_cell_alpha_beta_gamma_coordinates_and_panels(ranks,
                                          coarse_model,
                                          setup_coarse_cell_vertices_alpha_beta_coordinates(),
                                          collect(1:NPANELS),
                                          ptr_pXest_connectivity,
                                          ptr_new_pXest,
                                          ptr_pXest_ghost)

    GridapP4est.pXest_ghost_destroy(pXest_type,ptr_pXest_ghost)
    GridapP4est.pXest_lnodes_destroy(pXest_type,ptr_pXest_lnodes)

    pXest_refinement_rule_type = GridapP4est.PXestHorizontalRefinementRuleType()

    extruded_ref_coarsen_flags=
       map(partition(get_cell_gids(parametric_model.octree_dmodel.dmodel)),_refinement_and_coarsening_flags) do indices, flags
      similar(flags, length(local_to_global(indices)))
    end

    GridapP4est._extrude_refinement_and_coarsening_flags!(extruded_ref_coarsen_flags,
                                              _refinement_and_coarsening_flags,
                                              parametric_model.octree_dmodel.ptr_pXest,
                                              ptr_new_pXest)

    stride = GridapP4est.pXest_stride_among_children(pXest_type,
                                         pXest_refinement_rule_type,
                                         parametric_model.octree_dmodel.ptr_pXest)

    adaptivity_glue = GridapP4est._compute_fine_to_coarse_model_glue(pXest_type,
                                                       pXest_refinement_rule_type,
                                                       parametric_model.octree_dmodel.parts,
                                                       parametric_model.octree_dmodel.dmodel,
                                                       fmodel,
                                                       extruded_ref_coarsen_flags,
                                                       stride)

    adaptive_models = map(local_views(parametric_model.octree_dmodel),
                           local_views(fmodel),
                           adaptivity_glue) do model, fmodel, glue
           Gridap.Adaptivity.AdaptedDiscreteModel(fmodel,model,glue)
   end
   fmodel = GridapDistributed.GenericDistributedDiscreteModel(adaptive_models,get_cell_gids(fmodel))
   ref_model = OctreeDistributedDiscreteModel(3,3,
                                              ranks,
                                              fmodel,
                                              non_conforming_glue,
                                              coarse_model,
                                              ptr_pXest_connectivity,
                                              ptr_new_pXest,
                                              pXest_type,
                                              parametric_model.octree_dmodel.pXest_refinement_rule_type,
                                              false,
                                              parametric_model)

   # Build the proc-local ParametricDiscreteModels
   parametric_models = _setup_parametric_models(ref_model,
                                                cell_coordinates,
                                                panels)

   adaptive_models = map(parametric_models,
                         local_views(parametric_model.parametric_dmodel),
                         local_views(ref_model.dmodel)) do parametric_model,
                                                           parametric_model_parent,
                                                           octree_dmodel_adapted_model
        Gridap.Adaptivity.AdaptedDiscreteModel(parametric_model,
                                               parametric_model_parent,
                                               get_adaptivity_glue(octree_dmodel_adapted_model))
   end
   generic_dmodel =
      GenericDistributedDiscreteModel(adaptive_models, get_cell_gids(ref_model.dmodel))

   Parametric3DOctreeDistributedDiscreteModel(ref_model, generic_dmodel)
end
