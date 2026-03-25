"""
solve scalar laplacian in mixed form
u + ∇ᵧ(φ) = 0
∇ᵧ⋅u = f
where f = -Δφ
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

function petsc_mumps_setup(ksp)
  pc       = Ref{GridapPETSc.PETSC.PC}()
  mumpsmat = Ref{GridapPETSc.PETSC.Mat}()
  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPPREONLY)
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCLU)
  @check_error_code GridapPETSc.PETSC.PCFactorSetMatSolverType(pc[],GridapPETSc.PETSC.MATSOLVERMUMPS)
  @check_error_code GridapPETSc.PETSC.PCFactorSetUpMatSolverType(pc[])
  @check_error_code GridapPETSc.PETSC.PCFactorGetMatrix(pc[],mumpsmat)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  4, 1)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 28, 2)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 29, 1)
  # @check_error_code GridapPETSc.PETSC.MatMumpsSetCntl(mumpsmat[],  1, 0.00001)
  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
end


function f(p)
  function _f(α)
    xyz = ForwardMap(p)(α)
    θϕr   = xyz2θϕr(xyz)
    sin(θϕr[2])
  end
end



function launch_hodge_laplacian(ranks,n_ref,p_fe::Int,dir::String,return_vtk=1)

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("Hodge")

  octree3_model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
    num_horizontal_uniform_refinements=n_ref,
    num_vertical_uniform_refinements=n_ref);
  panel_model = octree3_model.parametric_dmodel

  GridapPETSc.Init()
  ls = PETScLinearSolver(petsc_mumps_setup)

  hodge_laplacian_scalar(panel_model,p_fe,dir,f,ls,Bool(return_vtk))

  GridapPETSc.Finalize()
  GridapPETSc.gridap_petsc_gc()

  i_am_main(ranks) && println("--DONE--")

end

function hodge_laplacian_scalar(
  panel_model::GridapDistributed.GenericDistributedDiscreteModel{3,3},
  p_fe::Int,dir::String,f::Function,ls=LUSolver(),return_vtk=false)

  ranks = get_ranks(panel_model)
  lvl_h = nref(nc_horizontal(panel_model))
  lvl_v = nref(nc_vertical(panel_model))
  i_am_main(ranks) && println("nref_h = $lvl_h; nref_v = $lvl_v; p_fe = $p_fe")

  degree = 5*(p_fe + 1)
  if p_fe == 0
    degree = 10
  end
  @check degree > 0 "Zero quad!!"

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)
  dΩ_error = Measure(Ω_error,2*degree)

  tags = ["top_boundary", "bottom_boundary"]
  Γ = BoundaryTriangulation(panel_model,tags=tags)
  dΓ = Measure(Γ,degree)
  nΓ = get_normal_vector(Γ)

  f_panel_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  sigma_cf = panelwise_cellfield(contr_gradf(f),Ω_panel,panel_ids)
  sdiv_cf =  panelwise_cellfield(surfdiv(contr_gradf(f)),Ω_panel,panel_ids)
  slap_panel_cf =  panelwise_cellfield(surflap(f),Ω_panel,panel_ids)
  rhs = -slap_panel_cf

  i_am_main(ranks) && println("Check sdiv(sgrad) against slap: ", l2(sdiv_cf-slap_panel_cf,dΩ) )

  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  # FE spaces
  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  biform1((u,p),(v,q)) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
  biform2((u,p),(v,q)) =  ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) )  )dΩ

  biformX((u,p),(v,q)) = biform1((u,p),(v,q)) + biform2((u,p),(v,q))
  liformX((v,q)) = ∫( (rhs*q)*meas_cf )dΩ + ∫( -f_panel_cf*(v⋅nΓ)*meas_cf )dΓ



  op = AffineFEOperator(biformX,liformX,X,Y)

  # A = get_matrix(op)
  # eigvals(Array(partition(A).item))
  # x = Gridap.Algebra.allocate_in_domain(A); fill!(x,1.0)
  # partition(A*x).item

  uh,ph = solve(ls,op)

  _e = f_panel_cf - ph
  el2_p = sqrt(sum(∫( (_e*_e)*meas_cf  )dΩ_error))

  _e = (covarient_basis_cf⋅uh) - (covarient_basis_cf⋅-sigma_cf) ### u = -∇p
  el2_u = sqrt(sum(∫( (_e⋅(metric_cf ⋅_e))*meas_cf  )dΩ_error))

  if return_vtk
    cellfields =  ["u"=>covarient_basis_cf⋅-sigma_cf,
    "uh"=>covarient_basis_cf⋅uh,
    "eu"=> (covarient_basis_cf⋅uh) - (covarient_basis_cf⋅-sigma_cf),
    "ph"=>ph, "p"=>f_panel_cf, "e"=>ph-f_panel_cf
                  ]
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl_h)_p$p_fe",
            cellfields=cellfields,
            append=false,geo_map= geo_map_func(Ω_panel))
  end

  ### convergence output for DrWatson
  dir_convergence = dir*"/convergence"
  (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

  n = nc(panel_model)
  n_h = nc_horizontal(panel_model)
  n_v = _nc_vertical(panel_model)
  dxx = dx(panel_model)
  dxH = dx_horizontal(panel_model)
  dxV = dx_vertical(panel_model)
  output = @strdict el2_u el2_p n n_h n_v dxx dxH dxV p_fe lvl_h lvl_v
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("hodge_laplacian_scalar_nrefh$(lvl_h)_nrefv$(lvl_v)_p$p_fe.jld2")), output)

  return el2_u, el2_p, false

end
