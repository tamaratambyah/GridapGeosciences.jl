"""
test mass conservation by solving
u + ∇ᵧ(φ) = u₀
∇ᵧ⋅u = 0
for arbitray vector field u₀
"""

using MPI
using PartitionedArrays

using DrWatson
using Gridap
using GridapDistributed
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra
using Gridap.ReferenceFEs, Gridap.Polynomials, Gridap.CellData

using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Test
using LinearAlgebra

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

function f(p)
  function _f(α)
    xyz = ForwardMap(p)(α)
    θϕr   = xyz2θϕr(xyz)
    sin(θϕr[2])
  end
end

p_fe = 1
ls = LUSolver()

#####   SERIAL 2D MODEL
# model0 = coarse_parametric_model()
# model = Gridap.Adaptivity.refine(model0)
# panel_model = model


##### DISTRIBUTED 2D MODEL
omodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=1)
panel_model = omodel.parametric_dmodel

##### DISTRIBUTED 3D MODEL
# octree3_model = Parametric3DOctreeDistributedDiscreteModel(ranks;
#         num_horizontal_uniform_refinements=1,
#         num_vertical_uniform_refinements=1);
# panel_model = octree3_model.parametric_dmodel


panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,2*(p_fe+1))

vec_contra_cf = panelwise_cellfield(contr_gradf(f),Ω_panel,panel_ids)


metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)
covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

# Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2,constraint=:zeromean)
Q0 = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
Q = Gridap.FESpaces.ZeroMeanFESpace(Q0,dΩ)
P = TrialFESpace(Q)

V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
U = TrialFESpace(V)

Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

biform1((u,p),(v,q)) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
biform2((u,p),(v,q)) =  ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) )  )dΩ

biformX((u,p),(v,q)) = biform1((u,p),(v,q)) + biform2((u,p),(v,q))
liformX((v,q)) = ∫( vec_contra_cf⋅(metric_cf⋅v)*meas_cf )dΩ

op = AffineFEOperator(biformX,liformX,X,Y)

A = get_matrix(op)
eigvals(Array(partition(A).item))



# uh,ph = solve(LUSolver(),op)

# # mass conservation errors
# s_div = sum(∫(  divergence(meas_cf*uh) )dΩ)
# s_div0 = sum(∫(  divergence(meas_cf*vec_contra_cf) )dΩ)
# panel_div = sum(∫(  divergence(uh) )dΩ)

# if return_vtk
#   lvl = nref(nc(panel_model))
#   cell_geo_map = geo_map_func(Ω_panel)

#   u_proj = covarient_basis_cf ⋅ vec_contra_cf
#   u_projh = covarient_basis_cf ⋅ uh

#   panel_cfs = [ u_proj, u_projh, ph]
#   labels = ["u0", "u_projh", "p"]

#   cellfields = map((x,y) -> x=>y, labels,panel_cfs)
#   writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",cellfields=cellfields,append=false,geo_map=cell_geo_map)
# end

# println("initial divergence: $s_div0")

# return abs(s_div), false,false
