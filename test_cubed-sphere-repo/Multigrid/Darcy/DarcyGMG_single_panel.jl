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
using GridapGeosciences


function get_patch_smoothers(sh,biform,qdegree)
  nlevs = num_levels(sh)
  smoothers = map(view(sh,1:nlevs-1)) do shl
    model = get_model(shl)
    ptopo = Geometry.PatchTopology(ReferenceFE{0},model)
    space = get_fe_space(shl)
    Ω  = Geometry.PatchTriangulation(model,ptopo)
    dΩ = Measure(Ω,qdegree)

    ap = (u,v) -> biform(u,v,dΩ)
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
  return (u,v) -> biform(u,v,dΩ)
end


MPI.Init()
np = MPI.Comm_size(MPI.COMM_WORLD)
parts = distribute_with_mpi(LinearIndices((np,)))
np_per_level = [(1,1),(1,1)]
ncell = (6,6)
Dc = length(ncell)
domain =  (-π/4,π/4,-π/4,π/4)

mh = CartesianModelHierarchy(parts,np_per_level,domain,ncell)
model = get_model(mh,1)

pid = 1
p_fe = 1
qdegree = 8*(p_fe+1)


h = panel_to_cartesian(h₀(0.0))
vX = panel_to_cartesian(tangent_vec(u₀(0.0)))

u_exact(x) = contra_v(vX)(1)(x)
p_exact(x) = h(1)(x)

# manufacture rhs functions
f(x) = u_exact(x) + inv_metric(pid)(x)⋅( ∇(p_exact)(x) )

tests_u  = TestFESpace(mh, ReferenceFE(raviart_thomas,Float64,p_fe);conformity=:HDiv,dirichlet_tags=["boundary"]);
trials_u = TrialFESpace(tests_u,[u_exact]);
U, V = get_fe_space(trials_u,1), get_fe_space(tests_u,1)
Q = TestFESpace(model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)

mfs = Gridap.MultiField.BlockMultiFieldStyle()
X = MultiFieldFESpace([U,Q];style=mfs)
Y = MultiFieldFESpace([V,Q];style=mfs)

α = 1.e2

metric_cf(x) = metric(pid)(x)
meas_cf(x) = sqrtg(pid)(x)
grad_meas_cf(x) = grad_meas(pid)(x)

inv_meas(x) = α*(sqrtg(pid)(x))^(-1)


graddiv(u,v,dΩ) = ∫( ( inv_meas)*( (v⋅grad_meas_cf + meas_cf*(∇⋅v) )*(u⋅grad_meas_cf + meas_cf*(∇⋅u) ) ))dΩ
biform_u(u,v,dΩ) = ∫( u⋅(metric_cf⋅v)*meas_cf )dΩ + graddiv(u,v,dΩ)
biform((u,p),(v,q),dΩ) = biform_u(u,v,dΩ) - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ - ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) ) )dΩ
liform((v,q),dΩ) = ∫( f⋅(metric_cf⋅v)*meas_cf )dΩ

Ω = Triangulation(model)
dΩ = Measure(Ω,qdegree)


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
  maxiter=10,mode=:preconditioner,verbose=1
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
function Gridap.CellData.get_triangulation(a::GridapDistributed.DistributedMultiFieldCellField)
  trians = map(get_triangulation,a.field_fe_fun)
  # @check all(map(t -> t === first(trians), trians))
  return first(trians)
end

xh = FEFunction(X,x)
uh,ph = xh

e = uh-u_exact
sqrt(sum(∫( e⋅e )dΩ )) # 1e-8

bb_cf = CellField(covariant_basis(1),Ω)
_e = bb_cf ⋅(e)
sqrt(sum(∫( ( _e⋅_e )*meas_cf)dΩ )) # 1e-8


using DrWatson
dir = datadir("Multigrid")
!isdir(dir) && mkdir(dir)
writevtk(Ω,dir*"/darcy_flat_alpha$α",cellfields= ["uh"=>bb_cf⋅uh, "ph"=>ph, "e"=>bb_cf⋅(uh-u_exact),
"u"=>bb_cf⋅u_exact,"divu"=>surfdiv(contra_v(vX))(1) ],append=false)
