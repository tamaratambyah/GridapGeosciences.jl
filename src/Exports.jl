macro publish(mod,name)
  quote
    using GridapGeosciences.$mod: $name; export $name
  end
end

@publish Adaptivity refine


@publish Fields ForwardMapPanel1
@publish Fields MatMultField
@publish Fields MyAffineField

@publish Geometry get_panel_ids
@publish Geometry geo_map_func
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

# @publish Helpers f_sin
# @publish Helpers fθϕ
# @publish Helpers f_XYZ
@publish Helpers rho
@publish Helpers rho3
@publish Helpers sqrtg
@publish Helpers _sqrtg
@publish Helpers detg
@publish Helpers grad_meas
@publish Helpers analytic_metric
@publish Helpers analytic_inv_metric
@publish Helpers _analytic_inv_metric

@publish Helpers analytic_perp_matrix
@publish Helpers surflap
@publish Helpers surfdiv
@publish Helpers contr_gradf
@publish Helpers sgrad

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

@publish Helpers forward_map
@publish Helpers forward_jacobian
@publish Helpers covarient_basis
@publish Helpers forward_pinv_jacobian
@publish Helpers contravariant_basis

@publish Distributed ParametricOctreeDistributedDiscreteModel
@publish Distributed DistributedParametricDiscreteModel
@publish Distributed DistributedAdaptedParametricDiscreteModel
@publish Distributed geo_map_func
@publish Distributed panelwise_cellfield

@publish Distributed writevtk
@publish Distributed createvtk
@publish Distributed write_vtk_file
@publish Distributed create_vtk_file
@publish Distributed create_pvtk_file
