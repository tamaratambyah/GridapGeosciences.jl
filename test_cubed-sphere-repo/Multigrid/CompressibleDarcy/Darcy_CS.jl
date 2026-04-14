"""
Solve the following on cubed sphere
u   + ∇p  = f
Cp  + ∇⋅u = g
where C = 1/c
"""

using Test
using LinearAlgebra
using FillArrays, BlockArrays

using Gridap
using Gridap.ReferenceFEs, Gridap.Algebra, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.MultiField, Gridap.Algebra
using PartitionedArrays
using GridapDistributed

using GridapSolvers
using GridapSolvers.LinearSolvers, GridapSolvers.MultilevelTools, GridapSolvers.PatchBasedSmoothers
using GridapSolvers.BlockSolvers: LinearSystemBlock, BiformBlock, BlockTriangularSolver
using MPI
using GridapP4est
using GridapGeosciences



function get_patch_smoothers(sh,biform,qdegree)
  nlevs = num_levels(sh)
  smoothers = map(view(sh,1:nlevs-1)) do shl
    model = get_model(shl)
    ptopo = Geometry.PatchTopology(ReferenceFE{0},model)
    space = get_fe_space(shl)
    Ω  = Geometry.PatchTriangulation(model,ptopo)
    dΩ = Measure(Ω,qdegree)

    panel_ids = get_panel_ids(model)
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

  panel_ids = get_panel_ids(model)
  metric_cf = panelwise_cellfield(metric,Ω,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω,panel_ids)
  grad_meas_cf = panelwise_cellfield(grad_meas,Ω,panel_ids)

  return (u,v) -> biform(u,v,dΩ,metric_cf,meas_cf,grad_meas_cf)
end

function Gridap.CellData.get_triangulation(a::GridapDistributed.DistributedMultiFieldCellField)
  trians = map(get_triangulation,a.field_fe_fun)
  # @check all(map(t -> t === first(trians), trians))
  return first(trians)
end


MPI.Init()
np = MPI.Comm_size(MPI.COMM_WORLD)
ranks = distribute_with_mpi(LinearIndices((np,)))


p_fe = 1
n_gmg_lvls = 2

model0 = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=1)
mh = ModelHierarchy(model0,n_gmg_lvls)

model = get_model(mh,1)
panel_ids = get_panel_ids(model)
Ω = Triangulation(model)
qdegree = 2*(p_fe+2)
dΩ = Measure(Ω,qdegree)

tests_u = TestFESpace(mh,ReferenceFE(raviart_thomas,Float64,p_fe);conformity=:Hdiv);
trials_u = TrialFESpace(tests_u);
U = get_fe_space(trials_u,1)
V = get_fe_space(tests_u,1)
Q = TestFESpace(model,ReferenceFE(lagrangian,Float64,p_fe);conformity=:L2)

mfs = Gridap.MultiField.BlockMultiFieldStyle()
X = MultiFieldFESpace([U,Q];style=mfs)
Y = MultiFieldFESpace([V,Q];style=mfs)

c = 0
α = 10

C = c ≠ 0 ? 1/c : 0.0

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
rhs_scalar = C*h_cf + sdiv_cf

metric_cf = panelwise_cellfield(metric,Ω,panel_ids)
meas_cf = panelwise_cellfield(sqrtg,Ω,panel_ids)
grad_meas_cf = panelwise_cellfield(grad_meas,Ω,panel_ids)


biform_u(u,v,dΩ,metric_cf,meas_cf,grad_meas_cf) = ( ∫( u⋅(metric_cf⋅v)*meas_cf )dΩ
                   + ∫( α*(1/meas_cf)*(v⋅grad_meas_cf + meas_cf*(∇⋅v) )*(u⋅grad_meas_cf + meas_cf*(∇⋅u) ) )dΩ
                  )
biform_p(p,q,dΩ,meas_cf) = ∫(C*(p*q)*meas_cf )dΩ
biform((u,p),(v,q),dΩ,metric_cf,meas_cf,grad_meas_cf) = (
                          biform_u(u,v,dΩ,metric_cf,meas_cf,grad_meas_cf)
                          + biform_p(p,q,dΩ,meas_cf)
                        - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
                        + ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) ) )dΩ
                          )
liform((v,q),dΩ,metric_cf,meas_cf,grad_meas_cf) = ∫( rhs_con_vector⋅(metric_cf⋅v)*meas_cf )dΩ + ∫( (q*rhs_scalar)*meas_cf  )dΩ

a(u,v) = biform(u,v,dΩ,metric_cf,meas_cf,grad_meas_cf)
l(v) = liform(v,dΩ,metric_cf,meas_cf,grad_meas_cf)
op = AffineFEOperator(a,l,X,Y)
A, b = get_matrix(op), get_vector(op);

#### solvers
biforms = map(mhl -> get_bilinear_form(mhl,biform_u,qdegree),mh)
smoothers = get_patch_smoothers(tests_u,biform_u,qdegree)
prolongations = setup_prolongation_operators(tests_u,qdegree;mode=:residual)
restrictions = setup_restriction_operators(
  tests_u,qdegree;mode=:residual,solver=CGSolver(JacobiLinearSolver())
)

gmg = GMGLinearSolver(
  trials_u,tests_u,biforms,
  prolongations,restrictions,
  pre_smoothers=smoothers,
  post_smoothers=smoothers,
  coarsest_solver=LUSolver(),
  maxiter=20,mode=:preconditioner,verbose=i_am_main(ranks),
  atol=1.0e-14, rtol=1.0e-08
)

##### solvers for the blocks of the preconditioner
solver_u = gmg
# solver_u = LUSolver()

solver_p = CGSolver(JacobiLinearSolver();maxiter=1000,atol=1e-14,rtol=1.e-8,verbose=1)
# solver_p = LUSolver()

#### preconditioner
bblocks  = [LinearSystemBlock() LinearSystemBlock();
            LinearSystemBlock() BiformBlock((p,q) -> ∫( (1.0/α + C)*(p*q)*meas_cf)dΩ,Q,Q)]
coeffs = [1.0 1.0;
          0.0 1.0]

P = BlockTriangularSolver(bblocks,[solver_u,solver_p],coeffs,:upper)
# P = JacobiLinearSolver()

##### Preconditioned external solver
ls = FGMRESSolver(20,P;maxiter=1000,atol=1e-14,rtol=1.e-8,verbose=true)
# ls = GMRESSolver(40;Pr=JacobiLinearSolver(),Pl=nothing,maxiter=2000,rtol=1.e-8,verbose=true)
# ls = LUSolver()
ns = numerical_setup(symbolic_setup(ls,A),A)

x = allocate_in_domain(A); fill!(x,0.0)
solve!(x,ns,b)



xh = FEFunction(X,x)
uh,ph = xh

dΩ_error = Measure(Ω,2*qdegree)
covariant_basis_cf = panelwise_cellfield(covariant_basis,Ω,panel_ids)
uh_proj = covariant_basis_cf⋅uh

eu = u_proj_cf - uh_proj
sqrt(sum(∫( (eu⋅eu)*meas_cf )dΩ_error )) # 1e-8

ep = ph-h_cf
sqrt(sum(∫( (ep⋅ep)*meas_cf )dΩ_error )) # 1e-5

using DrWatson
dir = datadir("Multigrid")
!isdir(dir) && mkdir(dir)
writevtk(Ω,dir*"/darcy_CS",
        cellfields= ["uh"=>uh_proj, "ph"=>ph, "eu"=>eu, "ep"=>ep],
        append=false,geo_map=latlon_geo_map_func(Ω))
