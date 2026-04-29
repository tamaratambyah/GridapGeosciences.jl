using MPI
using PartitionedArrays

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Test
using GridapPETSc

dir = datadir("dofs")
!isdir(dir) && mkdir(dir)


## pullback 3D vector to 3D chart
MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))


#################### sphere
n_ref_lvls = 0

omodel = GridapGeosciences.Distributed.CubedSphere2DParametricOctreeDistributedDiscreteModel(ranks;
  num_initial_uniform_refinements=n_ref_lvls)
_panel_model = omodel.parametric_dmodel
panel_model = _panel_model.models.item # get the serial model

p_fe = 1

panel_ids = get_panel_ids(panel_model)
Ω = Triangulation(panel_model)
dΩ = Measure(Ω,6)


function fV(p)
  function f(αβ)
    xyz = forward_map_2D(p)(αβ)
    VectorValue(-xyz[2],xyz[1],0.0)
  end
end

f_cf = ParametricCellField(contra_v(fV),Ω,panel_ids)

### True L2 values
RL2 = TestFESpace(panel_model, ReferenceFE(lagrangian,VectorValue{2,Float64},p_fe);conformity=:L2)
HL2 = TrialFESpace(RL2,f_cf)
f_h = interpolate(f_cf,HL2)
dL2 = collect(get_cell_dof_values(f_h))
map(x->round.(x,sigdigits=2),dL2)



# ### Incorrect H1 values
# RH1 = TestFESpace(panel_model, ReferenceFE(lagrangian,VectorValue{2,Float64},p_fe);conformity=:H1)
# HH1 = TrialFESpace(RH1,f_cf)
# f_h = interpolate(f_cf,HH1)
# dH1 = collect(get_cell_dof_values(f_h))
# map(x->round.(x,sigdigits=2),dH1)


R = TestFESpace(panel_model, ReferenceFE(lagrangian,VectorValue{2,Float64},p_fe);conformity=:L2)
dofs = get_cell_dof_ids(R)

## get the cell map for evaluation of the jacobians
coords =  QUAD.vertex_coords
cmap = get_cell_map(get_grid(panel_model))

### for each node, list:
###         the master cell
###         the vertex of the master (for evaluation of cmap)
###         the slave cell
###         the vertex of the slave (for evaluation of cmap)

masters = [6, 5, 6, 3, 6, 4, 5, 6]
m_vertex = [1,3,2,3,4,4,4,3]
slaves = [[1,5], [1,3], [1,2], [1,2], [2,4], [2,3], [3,4], [4,5]]
s_vertex = [ [1,1], [2,1], [3,1], [4,2], [3,2], [4,4], [2,3], [1,2]  ]

Mdofs = []
Coeffs = []
sDOF_to_dof = []

### get the matrix coeffients for each node, for each slave
for node = 1:8

  master_panel = masters[node]
  master_vertex = m_vertex[node]

  slave_panels = slaves[node]
  slave_vertexs = s_vertex[node]

  mat_coeffs = Matrix[]

  master_dof = dofs[master_panel][[master_vertex,master_vertex+4]]

  JM = forward_jacobian(master_panel)( cmap[master_panel](coords[master_vertex] ))

  ### For each slave, get the dofs and the mat coeffs
  for (slave_panel,slave_vertex) in zip(slave_panels,slave_vertexs)

    JSinv = forward_pinv_jacobian(slave_panel)(cmap[slave_panel](coords[slave_vertex] ))
    coeffs = JSinv⋅JM
    _mat_coeffs = Matrix(coeffs)
    push!(mat_coeffs,_mat_coeffs)
    _d = dofs[slave_panel][[slave_vertex,slave_vertex+4]]
    push!(sDOF_to_dof,_d)

  end
  master_dofs = fill(master_dof,4)
  push!(Mdofs,master_dofs)

  coeffs = [ mat_coeffs[1][1,:], mat_coeffs[1][2,:],
            mat_coeffs[2][1,:], mat_coeffs[2][2,:]]
  push!(Coeffs,coeffs)

end

### Now we have all constraints, build the FE space
sDOF_to_dofs = Gridap.Arrays.Table(vcat(Mdofs...))
sDOF_to_coeffs = Gridap.Arrays.Table(vcat(Coeffs...))
Rconstrained = Gridap.FESpaces.FESpaceWithLinearConstraints(vcat(sDOF_to_dof...), sDOF_to_dofs, sDOF_to_coeffs, R)

# Assemble the constrained system matrix
metric_cf = ParametricCellField(metric,Ω_panel,panel_ids)
meas_cf = ParametricCellField(sqrtg,Ω_panel,panel_ids)
a(u,v) = ∫( (u⋅(metric_cf⋅v))*meas_cf )dΩ
# Aconstrained=assemble_matrix(a,Rconstrained,Rconstrained)

l(v) = ∫( (vec_contra_cf⋅(metric_cf⋅v))*meas_cf )dΩ
op = AffineFEOperator(a,l,Rconstrained,Rconstrained)
fh = solve(LUSolver(),op)

_e = fh-vec_contra_h
e =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ))


# fh = interpolate(f_cf, Rconstrained)

dconstrained = collect(get_cell_dof_values(fh))

## Check that we obtain the L2 values
map(x->round.(x,sigdigits=2),dconstrained) .≈ map(x->round.(x,sigdigits=2),dL2)



# for (cell,i) in zip(slave_panels,slave_vertexs)
#   dof = [i,i+4]
#   _dcon = map(x->round.(x,sigdigits=2),dconstrained)[cell][dof]
#   _dL2 = map(x->round.(x,sigdigits=2),dL2)[cell][dof]
#  @test _dcon ≈ _dL2
# end
