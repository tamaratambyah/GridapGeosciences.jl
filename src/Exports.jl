macro publish(mod,name)
  quote
    using GridapGeosciences.$mod: $name; export $name
  end
end

@publish Adaptivity refine


@publish Fields MatMultField
@publish Fields MyAffineField
@publish Fields Cartesian2SphereicalMap
@publish Fields Cartesian2SphereicalMap3D

@publish Geometry get_panel_ids
@publish Geometry geo_map_func
@publish Geometry latlon_geo_map_func
@publish Geometry pullback_area_form
@publish Geometry pushforward_normal
@publish Geometry get_facet_normal
@publish Geometry get_mapped_facet_normal
@publish Geometry BoundaryTriangulation
@publish Geometry generate_ptr
@publish Geometry coarse_cube_model
@publish Geometry coarse_parametric_model
@publish Geometry R1p
@publish Geometry A_cube2panel
@publish Geometry A_panel2cube
@publish Geometry b_panel2cube
@publish Geometry ParametricDiscreteModel
@publish Geometry _pushforward_normal
@publish Geometry _pullback_area_form
@publish Geometry get_forward_map_generator
@publish Geometry get_radius
@publish Geometry get_thickness
@publish Geometry panelwise_cellfield

@publish ODEs DAEFEOperator

@publish Visualisation make_pvd
@publish Visualisation writevtk
@publish Visualisation createvtk
@publish Visualisation write_vtk_file
@publish Visualisation create_vtk_file
@publish Visualisation mapped_vtkpoints

@publish Helpers xyz2θϕr
@publish Helpers θϕ2xyz
@publish Helpers spherical_to_cartesian_matrix
@publish Helpers xyz2θϕ

@publish Helpers _sqrtg
@publish Helpers sqrtg
@publish Helpers detg
@publish Helpers grad_meas
@publish Helpers metric
@publish Helpers inv_metric


@publish Helpers perp_matrix
@publish Helpers surflap
@publish Helpers surfdiv
@publish Helpers contr_gradf
@publish Helpers sgrad
@publish Helpers g_star
@publish Helpers panel_to_cartesian
@publish Helpers panel_to_latlon

@publish Helpers vector_length
@publish Helpers normal_vec
@publish Helpers tangent_vec
@publish Helpers contra_v
@publish Helpers contra_v_comp
@publish Helpers contra_v_perp
@publish Helpers projection_v
@publish Helpers normal_vector_from_basis
@publish Helpers contra_v_comp3D
@publish Helpers contra_v_perp3D
@publish Helpers piola

@publish Helpers ForwardMap
@publish Helpers forward_jacobian
@publish Helpers covariant_basis
@publish Helpers forward_pinv_jacobian

@publish Helpers p_convergence_auto_test
@publish Helpers h_convergence_auto_test
@publish Helpers get_refined_models
@publish Helpers get_distributed_refined_models
@publish Helpers get_octree_refined_models
@publish Helpers get_3D_octree_refined_models
@publish Helpers nref
@publish Helpers nc
@publish Helpers nc_horizontal
@publish Helpers dx
@publish Helpers dx_horizontal
@publish Helpers convergence_rate

@publish Distributed ParametricOctreeDistributedDiscreteModel
@publish Distributed Parametric3DOctreeDistributedDiscreteModel
@publish Distributed DistributedParametricDiscreteModel
@publish Distributed geo_map_func
@publish Distributed latlon_geo_map_func
@publish Distributed panelwise_cellfield

@publish Distributed writevtk
@publish Distributed createvtk
@publish Distributed write_vtk_file
@publish Distributed create_vtk_file
@publish Distributed create_pvtk_file

@publish Distributed _make_pvd_distributed
@publish Distributed distributed_panel_ids
@publish Distributed DistributedAdaptivityGlue
@publish Distributed get_distributed_panel_model
@publish Distributed get_panel_ids
@publish Distributed get_owned_panel_ids
@publish Distributed get_skel_panel_ids
# @publish Distributed BoundaryTriangulation
@publish Distributed pullback_area_form
@publish Distributed pushforward_normal
@publish Distributed get_forward_map_generator
@publish Distributed get_radius
@publish Distributed get_thickness

@publish MultilevelTools ModelHierarchy
@publish MultilevelTools adapt_model
