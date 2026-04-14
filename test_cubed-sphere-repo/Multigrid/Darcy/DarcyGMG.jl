using Test
using LinearAlgebra
using FillArrays, BlockArrays

using Gridap
using Gridap.ReferenceFEs, Gridap.Algebra, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.MultiField, Gridap.Algebra
using PartitionedArrays
using GridapDistributed
using MPI
using GridapSolvers
using GridapSolvers.LinearSolvers, GridapSolvers.MultilevelTools, GridapSolvers.PatchBasedSmoothers
using GridapSolvers.BlockSolvers: LinearSystemBlock, BiformBlock, BlockTriangularSolver
using GridapGeosciences

function get_patch_smoothers(sh,biform,qdegree)
  nlevs = num_levels(sh)
  smoothers = map(view(sh,1:nlevs-1)) do shl
    model = get_model(shl)
    ptopo = Geometry.PatchTopology(ReferenceFE{0},model)
    space = get_fe_space(shl)
    Ω  = Geometry.PatchTriangulation(model,ptopo)
    dΩ = Measure(Ω,qdegree)

    panel_ids  = get_panel_ids(Ω)
    metric_cf = panelwise_cellfield(metric,Ω,panel_ids)
    meas_cf = panelwise_cellfield(sqrtg,Ω,panel_ids)
    grad_meas_cf = panelwise_cellfield(grad_meas,Ω,panel_ids)

    ap = (u,v) -> biform(u,v,dΩ,metric_cf,meas_cf,grad_meas_cf)
    solver = PatchBasedSmoothers.PatchSolver(
      ptopo, space, space, ap;
      assembly = :star,
      collect_factorizations = true,
      is_nonlinear = false
    )
    return RichardsonSmoother(solver,10,0.2)
  end
  return smoothers
end

function get_bilinear_form(mh_lev,biform,qdegree)
  model = get_model(mh_lev)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,qdegree)

  panel_ids  = get_panel_ids(Ω)
  metric_cf = panelwise_cellfield(metric,Ω,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω,panel_ids)
  grad_meas_cf = panelwise_cellfield(grad_meas,Ω,panel_ids)

  return (u,v) -> biform(u,v,dΩ,metric_cf,meas_cf,grad_meas_cf)
end


MPI.Init()
np = MPI.Comm_size(MPI.COMM_WORLD)
parts = distribute_with_mpi(LinearIndices((np,)))
np_per_level = [(1,1),(1,1)]
model0 = ParametricOctreeDistributedDiscreteModel(parts; num_initial_uniform_refinements=1)
n_ref_lvls = length(np_per_level)
mh = ModelHierarchy(model0,n_ref_lvls)

model = get_model(mh,1)



p_fe = 1
Ω = Triangulation(model)
qdegree = 6*(p_fe+1)
dΩ = Measure(Ω,qdegree)
panel_ids = get_panel_ids(model)


include("../../helpers.jl")
include("../../Geophysical/Williamson2Test.jl")
h = panel_to_cartesian(h₀(0.0))
vX = panel_to_cartesian(tangent_vec(u₀(0.0)))

h_cf = panelwise_cellfield(h,Ω,panel_ids)
u_proj_cf = panelwise_cellfield(projection_v(vX),Ω,panel_ids)
u_contra_cf = panelwise_cellfield(contra_v(vX),Ω,panel_ids)

covariant_basis_cf = panelwise_cellfield(covariant_basis,Ω,panel_ids)
sgrad_cf = panelwise_cellfield(sgrad(h),Ω,panel_ids)
pinvJ_cf = panelwise_cellfield(forward_pinv_jacobian,Ω,panel_ids)
sdiv_cf =  panelwise_cellfield(surfdiv(contra_v(vX)),Ω,panel_ids)

# manufacture rhs functions
rhs_vector = u_proj_cf + sgrad_cf
rhs_con_vector = pinvJ_cf ⋅ rhs_vector # exact contravariant component

tests_u  = TestFESpace(mh, ReferenceFE(raviart_thomas,Float64,p_fe);conformity=:HDiv);
trials_u = TrialFESpace(tests_u);
U, V = get_fe_space(trials_u,1), get_fe_space(tests_u,1)
Q =  TestFESpace(model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)

mfs = Gridap.MultiField.BlockMultiFieldStyle()
X = MultiFieldFESpace([U,Q];style=mfs)
Y = MultiFieldFESpace([V,Q];style=mfs)

α = 1.e2

metric_cf = panelwise_cellfield(metric,Ω,panel_ids)
meas_cf = panelwise_cellfield(sqrtg,Ω,panel_ids)
grad_meas_cf = panelwise_cellfield(grad_meas,Ω,panel_ids)

graddiv(u,v,dΩ,meas_cf,grad_meas_cf) = ∫( α*(1/meas_cf)*(v⋅grad_meas_cf + meas_cf*(∇⋅v) )*(u⋅grad_meas_cf + meas_cf*(∇⋅u) ) )dΩ
biform_u(u,v,dΩ,metric_cf,meas_cf,grad_meas_cf) = ∫( u⋅(metric_cf⋅v)*meas_cf )dΩ + graddiv(u,v,dΩ,meas_cf,grad_meas_cf)
biform((u,p),(v,q),dΩ) = biform_u(u,v,dΩ,metric_cf,meas_cf,grad_meas_cf) - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ - ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) ) )dΩ
liform((v,q),dΩ) = ∫( rhs_con_vector⋅(metric_cf⋅v)*meas_cf )dΩ

a(u,v) = biform(u,v,dΩ)
l(v) = liform(v,dΩ)
op = AffineFEOperator(a,l,X,Y)
A, b = get_matrix(op), get_vector(op);

biforms = map(mhl -> get_bilinear_form(mhl,biform_u,qdegree),mh)
smoothers = get_patch_smoothers(
  tests_u,biform_u,qdegree
)
prolongations = setup_prolongation_operators(
  tests_u,qdegree;mode=:residual
)
restrictions = setup_restriction_operators(
  tests_u,qdegree;mode=:residual,solver=CGSolver(JacobiLinearSolver())
)

gmg = GMGLinearSolver(
  trials_u,tests_u,biforms,
  prolongations,restrictions,
  pre_smoothers=smoothers,
  post_smoothers=smoothers,
  coarsest_solver=LUSolver(),
  maxiter=15,mode=:preconditioner,verbose=1
)

solver_u = gmg
solver_p = CGSolver(JacobiLinearSolver();maxiter=20,atol=1e-14,rtol=1.e-6,verbose=1)
solver_u.log.depth = 2
solver_p.log.depth = 4

bblocks  = [LinearSystemBlock() LinearSystemBlock();
            LinearSystemBlock() BiformBlock((p,q) -> ∫(-1.0/α*(p*q)*meas_cf)dΩ,Q,Q)]
coeffs = [1.0 1.0;
          0.0 1.0]
P = BlockTriangularSolver(bblocks,[solver_u,solver_p],coeffs,:upper)
solver = FGMRESSolver(20,P;atol=1e-14,rtol=1.e-14,verbose=1,maxiter=200)
ns = numerical_setup(symbolic_setup(solver,A),A)

x = allocate_in_domain(A); fill!(x,0.0)
solve!(x,ns,b)

r = allocate_in_range(A)
mul!(r,A,x)
r .-= b
@test norm(r) < 1.e-5
norm(r)

summary(solver_u.log)
summary(solver_p.log)
summary(solver.log)

gmg_iters = solver_u.log.num_iters
cg_iters = solver_p.log.num_iters
kylov_iters = solver.log.num_iters

function Gridap.CellData.get_triangulation(a::GridapDistributed.DistributedMultiFieldCellField)
  trians = map(get_triangulation,a.field_fe_fun)
  # @check all(map(t -> t === first(trians), trians))
  return first(trians)
end

xh = FEFunction(X,x)
uh,ph = xh

covariant_basis_cf = panelwise_cellfield(covariant_basis,Ω,panel_ids)
uh_proj = covariant_basis_cf⋅uh

e_u = sqrt(l2( (u_proj_cf - uh_proj),meas_cf,dΩ)) # error in physical velocity u
e_u = sqrt(l2( (uh - u_contra_cf),dΩ)) # error in physical velocity u

panel_cfs = [uh_proj, ph , u_proj_cf, u_proj_cf-covariant_basis_cf⋅uh, sdiv_cf]
labels = ["uh_proj", "ph", "u", "e", "sdivu"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)

cell_geo_map = geo_map_func(get_panel_ids(Ω))
latlon_cell_geo_map = latlon_geo_map_func(get_panel_ids(Ω))

dir = datadir("Multigrid")
!isdir(dir) && mkdir(dir)
writevtk(Ω,dir*"/darcy_sphere_alpha$(α)_ref",cellfields=cellfields,append=false,geo_map=latlon_cell_geo_map)
