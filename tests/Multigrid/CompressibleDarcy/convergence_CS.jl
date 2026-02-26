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

include("../../Geophysical/Williamson2Test.jl")
h = panel_to_cartesian(h₀(0.0))
vX = panel_to_cartesian(tangent_vec(u₀(0.0)))

### NOTE: n == n_ref_lvls
function convergence(ranks;c,α,n,order,iters,itu,itp,dir,return_vtk,simName)

  dir_convergence = dir*"/convergence"
  !isdir(dir_convergence) && mkdir(dir_convergence)

  # ensure no MPI task tries to generate the file before the main MPI task has
  # created the folder
  PartitionedArrays.barrier(ranks)

  C = c ≠ 0 ? 1/c : 0.0

  n_gmg_lvls = 1

  model0 = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=n)
  mh = ModelHierarchy(model0,n_gmg_lvls)

  model = get_model(mh,1)
  panel_ids = get_panel_ids(model)
  Ω = Triangulation(model)
  qdegree = 6*(order+2)
  dΩ = Measure(Ω,qdegree)

  tests_u = TestFESpace(mh,ReferenceFE(raviart_thomas,Float64,order);conformity=:Hdiv);
  trials_u = TrialFESpace(tests_u);
  U = get_fe_space(trials_u,1)
  V = get_fe_space(tests_u,1)
  Q = TestFESpace(model,ReferenceFE(lagrangian,Float64,order);conformity=:L2)

  mfs = Gridap.MultiField.BlockMultiFieldStyle()
  X = MultiFieldFESpace([U,Q];style=mfs)
  Y = MultiFieldFESpace([V,Q];style=mfs)

  ### panelwise cellfields
  h_cf = panelwise_cellfield(h,Ω,panel_ids)
  u_proj_cf = panelwise_cellfield(projection_v(vX),Ω,panel_ids)
  u_contra_cf = panelwise_cellfield(contra_v(vX),Ω,panel_ids)

  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω,panel_ids)
  sgrad_cf = panelwise_cellfield(sgrad(h),Ω,panel_ids)
  sgrad_contra_cf = panelwise_cellfield(contr_gradf(h),Ω,panel_ids)
  pinvJ_cf = panelwise_cellfield(forward_pinv_jacobian,Ω,panel_ids)
  sdiv_cf =  panelwise_cellfield(surfdiv(contra_v(vX)),Ω,panel_ids)

  # manufacture rhs functions
  # rhs_vector = u_proj_cf + sgrad_cf
  # rhs_con_vector = pinvJ_cf ⋅ rhs_vector # exact contravariant component
  rhs_con_vector = u_contra_cf + sgrad_contra_cf # exact contravariant component

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
  liform((v,q),dΩ,metric_cf,meas_cf,grad_meas_cf) = (
                   ∫( rhs_con_vector⋅(metric_cf⋅v)*meas_cf )dΩ
                 + ∫( (q*rhs_scalar)*meas_cf  )dΩ
                  )

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
    atol=1.0e-14, rtol=1.0e-12
  )

  cg = CGSolver(JacobiLinearSolver();maxiter=1000,atol=1e-14,rtol=1.e-12,verbose=i_am_main(ranks))

  ##### solvers for the blocks of the preconditioner
  solver_u = Bool(itu) ? gmg : LUSolver()
  solver_p = Bool(itp) ? cg  : LUSolver()

  #### preconditioner
  bblocks  = [LinearSystemBlock() LinearSystemBlock();
              LinearSystemBlock() BiformBlock((p,q) -> ∫( (1.0/α + C)*(p*q)*meas_cf)dΩ,Q,Q)]
  coeffs = [1.0 0.0;
            0.0 1.0]

  P = BlockTriangularSolver(bblocks,[solver_u,solver_p],coeffs,:upper)
  # P = JacobiLinearSolver()

  ##### Preconditioned external solver
  ls = FGMRESSolver(20,P;maxiter=iters,atol=1e-14,rtol=1.e-12,verbose=i_am_main(ranks))
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
  uh_proj = covarient_basis_cf⋅uh


  dΩ_error = Measure(Ω,2*qdegree)
  eu =  u_proj_cf - uh_proj
  Eu = sqrt(sum(∫( (eu⋅eu)*meas_cf )dΩ_error )) # 1e-8

  ep = ph-h_cf
  Ep = sqrt(sum(∫( (ep⋅ep)*meas_cf )dΩ_error )) # 1e-5

  if Bool(return_vtk)
    writevtk(Ω,dir*"/$(simName)",
      cellfields= ["uh"=>uh_proj, "ph"=>ph, "eu"=>eu, "ep"=>ep, "u"=>u_proj_cf],
      append=false,geo_map=latlon_geo_map_func(Ω))
  end

  output = @strdict Eu Ep gmg_iters cg_iters kylov_iters n order c α itu itp
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("$(simName).jld2")), output)

end
