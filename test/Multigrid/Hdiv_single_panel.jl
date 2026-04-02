using DrWatson
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

include("../convergence_tools.jl")
include("../Geophysical/Williamson2Test.jl")


# vX = panel_to_cartesian(tangent_vec(u₀(0.0)))
# u_exact(x) = contra_v(vX)(1)(x)
u_exact(x) = VectorValue(x[1],-x[2])
f(x) = u_exact(x)

tests  = TestFESpace(mh, ReferenceFE(raviart_thomas,Float64,p_fe);conformity=:HDiv,dirichlet_tags=["boundary"]);
trials = TrialFESpace(tests,u_exact);
spaces  = tests, trials

α = 1

metric_cf(x) = metric(pid)(x)
meas_cf(x) = sqrtg(pid)(x)
grad_meas_cf(x) = grad_meas(pid)(x)

inv_meas(x) = α*(sqrtg(pid)(x))^(-1)


graddiv(u,v,dΩ) = ∫( ( inv_meas)*( (v⋅grad_meas_cf + meas_cf*(∇⋅v) )*(u⋅grad_meas_cf + meas_cf*(∇⋅u) ) ))dΩ
biform_u(u,v,dΩ) = ∫( u⋅(metric_cf⋅v)*meas_cf )dΩ
biform(u,v,dΩ) = biform_u(u,v,dΩ) + graddiv(u,v,dΩ)
liform(v,dΩ) = ∫( f⋅(metric_cf⋅v)*meas_cf )dΩ

include("GMGTests.jl")

ctypes = [:v_cycle,:w_cycle,:f_cycle]
stypes = [:patch, :block]
ctype = ctypes[1]
stype = stypes[1]

t = 0
# gmg_driver_from_mats(t,parts,mh,spaces,biform,liform,u_exact,qdegree;ctype,stype)
# gmg_driver_from_weakform(t,parts,mh,spaces,biform,liform,u_exact,qdegree;ctype,stype)

  ctype=:v_cycle
  stype=:block
  ttype=:default

  # tests, trials = spaces
  restrictions, prolongations = get_transfers(ttype,mh,tests,biform,qdegree)
  smoothers = get_smoothers(stype,mh,tests,biform,qdegree)

  smatrices, A, b = compute_hierarchy_matrices(trials,tests,biform,liform,qdegree)
  gmg = GMGLinearSolver(
    smatrices,
    prolongations,restrictions,
    pre_smoothers=smoothers,
    post_smoothers=smoothers,
    coarsest_solver=LUSolver(),
    maxiter=1,mode=:preconditioner,
    cycle_type=ctype,
  )

  if ctype == :v_cycle
    solver = CGSolver(gmg;maxiter=20,atol=1e-14,rtol=1.e-6,verbose=i_am_main(parts))
  else
    solver = FGMRESSolver(5,gmg;maxiter=20,atol=1e-14,rtol=1.e-6,verbose=i_am_main(parts))
  end
  ns = numerical_setup(symbolic_setup(solver,A),A)

  # Solve

  x = pfill(0.0,partition(axes(A,2)))
  solve!(x,ns,b)


r = allocate_in_range(A)
mul!(r,A,x)
r .-= b
norm(r)

U    = get_fe_space(trials,1)
uh = FEFunction(U,x)

e = uh-u_exact

model = get_model(mh,1)
Ω = Triangulation(model)
dΩ = Measure(Ω,qdegree)

sqrt(sum(∫( e⋅e )dΩ )) # 1e-8

bb_cf = CellField(covarient_basis(1),Ω)
_e = bb_cf ⋅(e)
sqrt(sum(∫( ( _e⋅_e )*meas_cf)dΩ )) # 1e-8


using DrWatson
dir = datadir("Multigrid")
!isdir(dir) && mkdir(dir)
writevtk(Ω,dir*"/hdiv_flat_alpha$α",cellfields= ["uh"=>bb_cf⋅uh,  "e"=>bb_cf⋅(uh-u_exact),
"u"=>bb_cf⋅u_exact],append=false)
