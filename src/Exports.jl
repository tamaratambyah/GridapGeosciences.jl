macro publish(mod,name)
  quote
    using GridapGeosciences.$mod: $name; export $name
  end
end

@publish Adaptivity refine


@publish Fields MatMultField
@publish Fields MyAffineField
@publish Fields Cartesian2SphericalMap
@publish Fields ForwardMap
@publish Fields normal_vec

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
@publish Geometry _pushforward_normal
@publish Geometry _pullback_area_form
@publish Geometry get_forward_map_generator
@publish Geometry get_radius
@publish Geometry get_thickness
@publish Geometry ParametricCellField
@publish Geometry CubedSphereParametricDiscreteModel
@publish Geometry CubedSphere2DParametricDiscreteModel
@publish Geometry CubedSphere3DParametricDiscreteModel
@publish Geometry get_refined_models

@publish ODEs DAEFEOperator

@publish Visualisation make_pvd
@publish Visualisation writevtk_with_cell_geomap
@publish Visualisation createvtk_with_cell_geomap
@publish Visualisation write_vtk_file_with_cell_geomap
@publish Visualisation create_vtk_file_with_cell_geomap
@publish Visualisation mapped_vtkpoints

@publish Helpers xyz2θϕr
@publish Helpers θϕ2xyz
@publish Helpers xyz2θϕ

@publish Helpers sqrtg
@publish Helpers detg
@publish Helpers grad_meas
@publish Helpers metric
@publish Helpers inv_metric

@publish Helpers surflap
@publish Helpers surfdiv
@publish Helpers contr_gradf
@publish Helpers sgrad
@publish Helpers panel_to_cartesian

@publish Helpers tangent_vec
@publish Helpers contra_v
@publish Helpers piola

@publish Helpers forward_jacobian
@publish Helpers covariant_basis
@publish Helpers forward_pinv_jacobian

@publish Distributed CubedSphere2DParametricOctreeDistributedDiscreteModel
@publish Distributed CubedSphere3DParametricOctreeDistributedDiscreteModel
@publish Distributed CubedSphereParametricDistributedDiscreteModel
@publish Distributed CubedSphere2DParametricDistributedDiscreteModel
@publish Distributed CubedSphere3DParametricDistributedDiscreteModel
@publish Distributed geo_map_func
@publish Distributed latlon_geo_map_func
@publish Distributed ParametricCellField

@publish Distributed writevtk_with_cell_geomap
@publish Distributed createvtk_with_cell_geomap
@publish Distributed write_vtk_file_with_cell_geomap
@publish Distributed create_vtk_file_with_cell_geomap
@publish Distributed create_pvtk_file_with_cell_geomap

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

@publish ConvergenceTools p_convergence_auto_test
@publish ConvergenceTools h_convergence_auto_test
@publish ConvergenceTools get_distributed_refined_models
@publish ConvergenceTools get_octree_refined_models
@publish ConvergenceTools get_3D_octree_refined_models
@publish ConvergenceTools nref
@publish ConvergenceTools nc
@publish ConvergenceTools nc_horizontal
@publish ConvergenceTools dx
@publish ConvergenceTools dx_horizontal
@publish ConvergenceTools convergence_rate
