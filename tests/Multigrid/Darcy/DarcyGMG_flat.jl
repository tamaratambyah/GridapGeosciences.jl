
""" DarcyGMGApplication
- Flat space
- from GridapSolvers, develop branch
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
nc = (4,4)
Dc = length(nc)
domain =  (0,1,0,1)

mh = CartesianModelHierarchy(parts,np_per_level,domain,nc)
model = get_model(mh,1)

order = 2
qdegree = 2*(order+1)
reffe_u = ReferenceFE(raviart_thomas,Float64,order-1)
reffe_p = ReferenceFE(lagrangian,Float64,order-1;space=:P)

u_exact(x) = VectorValue(x[1]+x[2],-x[2])
p_exact(x) = 2.0*x[1]-1.0

tests_u  = TestFESpace(mh,reffe_u,dirichlet_tags=["boundary"]);
trials_u = TrialFESpace(tests_u,[u_exact]);
U, V = get_fe_space(trials_u,1), get_fe_space(tests_u,1)
Q = TestFESpace(model,reffe_p;conformity=:L2)

mfs = Gridap.MultiField.BlockMultiFieldStyle()
X = MultiFieldFESpace([U,Q];style=mfs)
Y = MultiFieldFESpace([V,Q];style=mfs)

α = 1.e2
f(x) = u_exact(x) + ∇(p_exact)(x)
graddiv(u,v,dΩ) = ∫(α*divergence(u)⋅divergence(v))dΩ
biform_u(u,v,dΩ) = ∫(v⊙u)dΩ + graddiv(u,v,dΩ)
biform((u,p),(v,q),dΩ) = biform_u(u,v,dΩ) - ∫(divergence(v)*p)dΩ - ∫(divergence(u)*q)dΩ
liform((v,q),dΩ) = ∫(v⋅f)dΩ

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
  maxiter=3,mode=:preconditioner,verbose=i_am_main(parts)
)

solver_u = gmg
solver_p = CGSolver(JacobiLinearSolver();maxiter=20,atol=1e-14,rtol=1.e-6,verbose=i_am_main(parts))
solver_p.log.depth = 2

bblocks  = [LinearSystemBlock() LinearSystemBlock();
            LinearSystemBlock() BiformBlock((p,q) -> ∫(-1.0/α*p*q)dΩ,Q,Q)]
coeffs = [1.0 1.0;
          0.0 1.0]
P = BlockTriangularSolver(bblocks,[solver_u,solver_p],coeffs,:upper)
solver = FGMRESSolver(20,P;atol=1e-14,rtol=1.e-10,verbose=i_am_main(parts))
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

using DrWatson
dir = datadir("Multigrid")
!isdir(dir) && mkdir(dir)
writevtk(Ω,dir*"/darcy_flat",cellfields= ["uh"=>uh, "ph"=>ph, "e"=>uh-u_exact],append=false)

e = uh-u_exact
sqrt(sum(∫( e⋅e )dΩ )) # 1e-8
