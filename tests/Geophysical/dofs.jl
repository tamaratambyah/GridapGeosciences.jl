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

omodel = GridapGeosciences.Distributed.ParametricOctreeDistributedDiscreteModel(ranks;
  num_initial_uniform_refinements=n_ref_lvls)
panel_model = omodel.parametric_dmodel

p_fe = 1

panel_ids = get_panel_ids(panel_model)
Ω = Triangulation(panel_model)
dΩ = Measure(Ω,6)

metric_cf = panelwise_cellfield(metric,Ω,panel_ids)
meas_cf = panelwise_cellfield(sqrtg,Ω,panel_ids)
covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω,panel_ids)

function fV(p)
  function f(αβ)
    xyz = forward_map_2D(p)(αβ)
    VectorValue(-xyz[2],0.0,0.0)
  end
end

f_cf = panelwise_cellfield(contra_v(fV),Ω,panel_ids)

### True L2 values
RL2 = TestFESpace(panel_model, ReferenceFE(lagrangian,VectorValue{2,Float64},p_fe);conformity=:L2)
HL2 = TrialFESpace(RL2,f_cf)
f_h = interpolate(f_cf,HL2)
dL2 = collect(get_cell_dof_values(f_h.fields.item_ref[]))
map(x->round.(x,sigdigits=2),dL2)

### Incorrect H1 values
RH1 = TestFESpace(panel_model, ReferenceFE(lagrangian,VectorValue{2,Float64},p_fe);conformity=:H1)
HH1 = TrialFESpace(RH1,f_cf)
f_h = interpolate(f_cf,HH1)
dH1 = collect(get_cell_dof_values(f_h.fields.item_ref[]))
map(x->round.(x,sigdigits=2),dH1)
