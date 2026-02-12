"""
Solve the following on periodic meshes
u   + ∇p  = f
Cp  + ∇⋅u = g
where C = 1/c
"""

using DataStructures
using DrWatson
using DataFrames

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


function adapt_model(ranks,model::OctreeDistributedDiscreteModel)
  cell_partition=get_cell_gids(model.dmodel)
  ref_flags=map(ranks,partition(cell_partition)) do rank,indices
      flags=zeros(Cint,length(indices))
      flags.=refine_flag
  end
  ref_model, adaptivity_glue = Gridap.Adaptivity.adapt(model,ref_flags)
  return ref_model, adaptivity_glue
end

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

function Gridap.CellData.get_triangulation(a::GridapDistributed.DistributedMultiFieldCellField)
  trians = map(get_triangulation,a.field_fe_fun)
  # @check all(map(t -> t === first(trians), trians))
  return first(trians)
end

u_exact(x) = VectorValue(x[2]*(1-x[2]),x[1]*(1-x[1]))
p_exact(x) = 1 + x[1]*(1-x[1])

# u_exact(x) = VectorValue(sin(2*π*x[2]),0.0)
# p_exact(x) = 1 + 0.1*sin(2*π*x[1])



function convergence(ranks;c,α,n,order,iters,itu,itp,dir,return_vtk,simName)

  dir_convergence = dir*"/convergence"
  !isdir(dir_convergence) && mkdir(dir_convergence)

  C = c ≠ 0 ? 1/c : 0.0

  coarse_model = CartesianDiscreteModel((0,1,0,1),(n,n),isperiodic=(true,true))
  dmodel = OctreeDistributedDiscreteModel(ranks, coarse_model)
  fmodel,glue = adapt_model(ranks,dmodel)

  mh = ModelHierarchy([fmodel.dmodel,dmodel.dmodel])
  model = get_model(mh,1)
  Ω = Triangulation(model)
  qdegree = 2*(order+2)
  dΩ = Measure(Ω,qdegree)

  tests_u = TestFESpace(mh,ReferenceFE(raviart_thomas,Float64,order);conformity=:Hdiv);
  trials_u = TrialFESpace(tests_u);
  U = get_fe_space(trials_u,1)
  V = get_fe_space(tests_u,1)
  Q = TestFESpace(model,ReferenceFE(lagrangian,Float64,order);conformity=:L2)

  mfs = Gridap.MultiField.BlockMultiFieldStyle()
  X = MultiFieldFESpace([U,Q];style=mfs)
  Y = MultiFieldFESpace([V,Q];style=mfs)


  f(x) = u_exact(x) + ∇(p_exact)(x)
  g(x) = C*p_exact(x) + (∇⋅u_exact)(x)
  biform_u(u,v,dΩ) = ∫(u⋅v)dΩ + ∫(α*(divergence(u)*divergence(v)) )dΩ
  biform_p(p,q,dΩ) = ∫(C*p*q )dΩ
  biform((u,p),(v,q),dΩ) = ( biform_u(u,v,dΩ) + biform_p(p,q,dΩ)
                            - ∫(divergence(v)*p)dΩ
                            + ∫(divergence(u)*q)dΩ
                            )
  liform((v,q),dΩ) = ∫(v⋅f)dΩ + ∫(q*g)dΩ

  a(u,v) = biform(u,v,dΩ)
  l(v) = liform(v,dΩ)
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

  cg = CGSolver(JacobiLinearSolver();maxiter=1000,atol=1e-14,rtol=1.e-8,verbose=i_am_main(ranks))

  ##### solvers for the blocks of the preconditioner
  solver_u = Bool(itu) ? gmg : LUSolver()
  solver_p = Bool(itp) ? cg  : LUSolver()

  #### preconditioner
  bblocks  = [LinearSystemBlock() LinearSystemBlock();
              LinearSystemBlock() BiformBlock((p,q) -> ∫( (1.0/α + C)*p*q)dΩ,Q,Q)]
  coeffs = [1.0 0.0;
            0.0 1.0]

  P = BlockTriangularSolver(bblocks,[solver_u,solver_p],coeffs,:upper)
  # P = JacobiLinearSolver()

  ##### Preconditioned external solver
  ls = FGMRESSolver(20,P;maxiter=iters,atol=1e-14,rtol=1.e-8,verbose=i_am_main(ranks))
  # ls = GMRESSolver(40;Pr=JacobiLinearSolver(),Pl=nothing,maxiter=2000,rtol=1.e-8,verbose=i_am_main(ranks))
  # ls = LUSolver()
  ns = numerical_setup(symbolic_setup(ls,A),A)

  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)

  gmg_iters = Bool(itu) ? solver_u.log.num_iters : 0
  cg_iters = Bool(itp) ? solver_p.log.num_iters : 0
  kylov_iters = ls.log.num_iters

  xh = FEFunction(X,x)
  uh,ph = xh

  dΩ_error = Measure(Ω,2*qdegree)
  eu = uh-u_exact
  Eu = sqrt(sum(∫( eu⋅eu )dΩ_error )) # 1e-8

  ep = ph-p_exact
  Ep = sqrt(sum(∫( ep⋅ep )dΩ_error )) # O(1)

  if Bool(return_vtk)
    writevtk(Ω,dir*"/$(simName)",
      cellfields= ["uh"=>uh, "ph"=>ph, "eu"=>eu, "ep"=>ep, "u"=>u_exact],
      append=false)
  end

  output = @strdict Eu Ep gmg_iters cg_iters kylov_iters n order c α itu itp
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("$(simName).jld2")), output)

end
