using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed
using DrWatson
using LinearAlgebra
using FillArrays

include("../convergence_tools.jl")


function interpolation(dpanel_model::GridapDistributed.GenericDistributedDiscreteModel{2,2},
                        p_fe::Int,dir::String,vecX::Function,return_vtk)

  panel_model = dpanel_model.models.item_ref[]

  Dc = num_cell_dims(panel_model)

  lvl = nref(nc(panel_model))
  println("p_fe = $(p_fe); nref = $lvl; Dc = $Dc")

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,4*p_fe)
  dΩ_error = Measure(Ω_panel,8*p_fe)
  panel_ids = get_panel_ids(panel_model)

  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  vec_contra_cf = panelwise_cellfield(contra_v(vecX),Ω_panel,panel_ids)
  vec_proj_cf = covarient_basis_cf⋅vec_contra_cf

  ### BEGING ALBERTO'S IMPLEMENTATION
  value_type=VectorValue{2,Float64}

  reffe = ReferenceFE(lagrangian,value_type,p_fe)
  lag_reffe = Gridap.ReferenceFEs.LagrangianRefFE(value_type,QUAD,p_fe)

  cell_map = get_cell_map(panel_model)
  reffe_node_coordinates = Gridap.ReferenceFEs.get_node_coordinates(lag_reffe)
  cell_wise_reffe_node_coordinates = lazy_map(evaluate, cell_map, Fill(reffe_node_coordinates, num_cells(panel_model)))

  V = FESpace(panel_model, reffe)

  panel_ids = get_panel_ids(panel_model)
  grid_topology = Gridap.Geometry.get_grid_topology(panel_model)
  grid = get_grid(panel_model)

  ndofs = Gridap.ReferenceFEs.num_dofs(lag_reffe)

  face_own_dofs  = Gridap.ReferenceFEs.get_face_own_dofs(lag_reffe)
  face_own_nodes = Gridap.ReferenceFEs.get_face_own_nodes(lag_reffe)

  S_matrices = Vector{Matrix{Float64}}(undef,num_cells(panel_model))

  Dc = num_dims(panel_model)

  # Loop over cells
  for cell=1:num_cells(panel_model)
      # Initialize transformation matrix S for the current cell
      S_matrices[cell] = Array{Float64}(I,(ndofs,ndofs))
      S = S_matrices[cell]

      current_cell_panel = panel_ids[cell]


      for i=0:Dc-1
          cell_faces = Gridap.Geometry.get_faces(grid_topology, Dc, i)

          # Warning: we must use the global IDs of the cells instead
          # of the proc-local IDs as we are doing now. With one processor
          # is ok, as both coincide.
          face_cells  = Gridap.Geometry.get_faces(grid_topology, i, Dc)

          num_faces_dim_i = num_faces(QUAD, i)
          offset = Gridap.ReferenceFEs.get_offset(QUAD,i)


          for j=1:num_faces_dim_i
              face_gid = cell_faces[cell][j]
              max_cell_gid = -1
              for cell_around in face_cells[face_gid]
                  if (cell_around>max_cell_gid)
                      max_cell_gid = cell_around
                  end
              end
              master_cell_panel = panel_ids[max_cell_gid]
              dof_lids_slave = face_own_dofs[offset+j]
              node_lids_slave = face_own_nodes[offset+j]
              if (current_cell_panel!=master_cell_panel && length(dof_lids_slave)>0)
                   # Vertex with local ID i within the current cell is a slave vertex.
                   # Determine local ID j within the master cell
                   pos_master=findfirst(x->x==face_gid,cell_faces[max_cell_gid])
                   node_lids_master = face_own_nodes[offset+pos_master]

                   for (inode, node_slave) in enumerate(node_lids_slave)

                       node_master = node_lids_master[inode]

                       JM = forward_jacobian(master_cell_panel)(cell_wise_reffe_node_coordinates[max_cell_gid][node_master])
                       JSinv = forward_pinv_jacobian(current_cell_panel)(cell_wise_reffe_node_coordinates[cell][node_slave])

                       #println("Cell $cell, vertex $i ($(cell_wise_reffe_node_coordinates[cell][i])) is slave. Master cell: $max_cell_gid vertex $j ($(cell_wise_reffe_node_coordinates[max_cell_gid][j]))")

                       coeffs = JSinv⋅JM

                       #println("coeffs=$coeffs")

                       dof_lids_slave_current_node = findall(x->x==node_slave,lag_reffe.reffe.dofs.dof_to_node)

                       # To-think: should we use here coeffs or its transpose?
                       S[dof_lids_slave_current_node,dof_lids_slave_current_node] .= Array(coeffs)
                   end
              end
          end
      end
      #println("S matrix for cell $cell (panel $current_cell_panel):")
      #display(S)
  end

  test_basis = get_fe_basis(V)
  test_basis_field_array = Gridap.CellData.get_data(test_basis)

  # VERY IMPORTANT: linear_combination works s.t. the i-th field in the output basis is the linear combination of the fields in input basis using the coefficients in the i-th COLUMN of the matrix. Thus, if we use S_matrices, we are actually building \phi = S^T \psi.
  transformed_test_basis_field_array =
      lazy_map(Gridap.Fields.linear_combination, S_matrices, test_basis_field_array)

  transformed_trial_basis_field_array=
     lazy_map(transpose, lazy_map(Gridap.Fields.linear_combination, S_matrices, test_basis_field_array))

  test_basis_transformed = Gridap.FESpaces.SingleFieldFEBasis(
      transformed_test_basis_field_array,
      Triangulation(panel_model), Gridap.FESpaces.TestBasis(), Gridap.CellData.ReferenceDomain())

  trial_basis_transformed = Gridap.FESpaces.SingleFieldFEBasis(
      transformed_trial_basis_field_array,
      Triangulation(panel_model), Gridap.FESpaces.TrialBasis(), Gridap.CellData.ReferenceDomain())

  a(u,v) = ∫( (u⋅(metric_cf⋅v))*meas_cf )dΩ
  l(v) = ∫( (vec_contra_cf⋅(metric_cf⋅v))*meas_cf )dΩ

  V = TestFESpace(panel_model, reffe; conformity=:H1)
  U = TrialFESpace(V)

  du=get_trial_fe_basis(U)
  dv=get_fe_basis(V)

  AK_dc_no_transformed = a(du,dv)
  AK_no_transformed = get_array(AK_dc_no_transformed)
  bK_dc_no_transformed = l(dv)
  bK_no_transformed = get_array(bK_dc_no_transformed)

  AK_dc_transformed = a(trial_basis_transformed, test_basis_transformed)
  AK_transformed = get_array(AK_dc_transformed)
  bK_dc_transformed = l(test_basis_transformed)
  bK_transformed = get_array(bK_dc_transformed)

  S_matrices_transposed = lazy_map(transpose, S_matrices)
  for i=1:num_cells(panel_model)
    @assert norm(S_matrices_transposed[i]*bK_no_transformed[i]-bK_transformed[i]) <1.0e-12
    @assert norm(S_matrices_transposed[i]*AK_no_transformed[i]*S_matrices[i]-AK_transformed[i]) < 1.0e-12
  end

  A=assemble_matrix(AK_dc_transformed, U, V)
  b=assemble_vector(bK_dc_transformed, V)
  udofs = A\b

  dirichlet_values = Float64[]

  cell_dofs = Gridap.FESpaces.scatter_free_and_dirichlet_values(V, udofs, dirichlet_values)

  fe_function_cell_field_array = lazy_map(Gridap.Fields.linear_combination,
                                          cell_dofs,
                                          lazy_map(Gridap.Fields.linear_combination, S_matrices, test_basis_field_array))
  fe_function_cell_field =
       Gridap.CellData.GenericCellField(fe_function_cell_field_array,Ω_panel,Gridap.CellData.ReferenceDomain())

  vec_contra_h = Gridap.FESpaces.SingleFieldFEFunction(fe_function_cell_field,cell_dofs,udofs,dirichlet_values,V)
  vec_l2proj_h = covarient_basis_cf ⋅vec_contra_h

  _e = vec_contra_cf - vec_contra_h
  el2_proj =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ_error))

  # Interpolation instead of L2 projection ...
  test_h = interpolate(vec_contra_cf, V)
  test_h_dofs = get_free_dof_values(test_h)
  dirichlet_values = Float64[]
  cell_dofs = Gridap.FESpaces.scatter_free_and_dirichlet_values(V, test_h_dofs, dirichlet_values)

  fe_function_cell_field_array = lazy_map(Gridap.Fields.linear_combination,
                                          cell_dofs,
                                          lazy_map(Gridap.Fields.linear_combination, S_matrices, test_basis_field_array))
  fe_function_cell_field =
       Gridap.CellData.GenericCellField(fe_function_cell_field_array,Ω_panel,Gridap.CellData.ReferenceDomain())

  vec_contra_h = Gridap.FESpaces.SingleFieldFEFunction(fe_function_cell_field,cell_dofs,test_h_dofs,dirichlet_values,V)
  vec_interp_h = covarient_basis_cf ⋅vec_contra_h

  test_h = interpolate(vec_contra_cf, V)
  test_h_dofs = get_free_dof_values(test_h)

  _e = vec_contra_cf - vec_contra_h
  el2_interp =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ_error))

  println("Error interp: ", el2_interp)
  println("Error proj: ", el2_proj)

  if return_vtk
    panel_cfs = [vec_proj_cf, vec_l2proj_h, vec_proj_cf-vec_l2proj_h,
                vec_interp_h, vec_interp_h-vec_proj_cf]
    labels = ["u_proj", "u_projh", "eproj",
              "u_int", "e_int"]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$(p_fe)",cellfields=cellfields,
          append=false,geo_map=latlon_geo_map_func(Ω_panel))
  end

  return  el2_interp,el2_proj,false

end

### must be in the tangent space of the sphere
vX(xyz) = VectorValue(-xyz[2], xyz[1], 0)
vecX = panel_to_cartesian(tangent_vec(vX))



MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))


n_ref_lvls = 4
ps = [2]

dir = datadir("InterpolationConvergence")
!isdir(dir) && mkdir(dir)

Dc = 2
models  = get_octree_refined_models(ranks,n_ref_lvls)


_dir = dir*"/vector_func_$(Dc)D_H1"
!isdir(_dir) && mkdir(_dir)
p_convergence_test(ranks,ps,models,interpolation,_dir,vecX,true)
plot_convergence_from_saved(_dir,"convergence",["Interp","L2Proj", ])
