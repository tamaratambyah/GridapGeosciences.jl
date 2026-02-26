
using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed
using LinearAlgebra
using FillArrays

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

parametric_octree_dmodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=0)

order=3
value_type=VectorValue{2,Float64}

reffe = ReferenceFE(lagrangian,value_type,order)
lag_reffe = Gridap.ReferenceFEs.LagrangianRefFE(value_type,QUAD,order)
panel_model = parametric_octree_dmodel.parametric_dmodel.models.item_ref[]

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

### must be in the tangent space of the sphere
vX(xyz) = VectorValue(-xyz[2], xyz[1], 0)
vecX = panel_to_cartesian(tangent_vec(vX))

Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,6)
panel_ids = get_panel_ids(panel_model)

metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

vec_contra_cf = panelwise_cellfield(contra_v(vecX),Ω_panel,panel_ids)
vec_proj_cf = covarient_basis_cf⋅vec_contra_cf

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

_e = vec_contra_cf - vec_contra_h
e =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ))

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


test_h = interpolate(vec_contra_cf, V)
test_h_dofs = get_free_dof_values(test_h)

# BEGIN: THIS PART OF CODE ONLY MAKES SENSE WITH SIX CELLS
# Master cell 3, vertex 3 
test_h_dofs[7]
test_h_dofs[8]

# Transforming to slave cell 1, vertex 4
S_matrices[1][[4,8],[4,8]]*test_h_dofs[[7,8]]

evaluate(test_basis_field_array[1][4] + 0.5*test_basis_field_array[1][8], Point(-0.5,0.5))
kk=lazy_map(Gridap.Fields.linear_combination, S_matrices, test_basis_field_array)
evaluate(kk[1][4] - 0.5*kk[1][8], Point(-0.5,0.5))
# END: THIS PART OF CODE ONLY MAKES SENSE WITH SIX CELLS


_e = vec_contra_cf - vec_contra_h
e =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ))