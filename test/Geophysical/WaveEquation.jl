"""
solve the linearised wave equation in steady form using manufactured solutions
u + ∇ᵧ(φ) = f₁
φ + ∇ᵧ⋅u = f₁
"""

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using GridapGeosciences
using Test

using GridapPETSc
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

include("Williamson2Test.jl")


function wave_solver(panel_model,
  p_fe::Int,dir::String,h::Function,vX::Function,ls=LUSolver(),return_vtk=false)

  ranks = get_ranks(panel_model)
  Dc = num_cell_dims(panel_model)
  lvl = nref(panel_model)

  i_am_main(ranks) && println("nref = $lvl; p_fe = $p_fe; Dc = $Dc")

  degree = 5*(p_fe+1)
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)
  dΩ_error = Measure(Ω_panel,2*degree)

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  # metric information
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
  u_cf = panelwise_cellfield(piola(vX),Ω_panel,panel_ids)
  u_proj_cf = covarient_basis_cf ⋅(1/meas_cf * u_cf  )

  p_int = interpolate(h_cf,P)
  u_int = interpolate(u_cf,U)

  biform1((u,p),(v,q)) = ∫( (u⋅ (metric_cf⋅v))*(1/meas_cf) )dΩ - ∫( p*(∇⋅v) )dΩ
  biform2((u,p),(v,q)) = ∫( (p*q)*meas_cf )dΩ + ∫( q*(∇⋅u) )dΩ
  biformX((u,p),(v,q)) = biform1((u,p),(v,q)) + biform2((u,p),(v,q))

  # manufacture rhs functions
  function get_liform(Dc::Int)

    # the manufactured solution is exactly the LHS operator
    _liformX((v,q)) = (
      ∫( (u_int⋅ (metric_cf⋅v))*(1/meas_cf) )dΩ
    + ∫( gradient(p_int)⋅v )dΩ # assume regularity to IBP
    + ∫( (p_int*q)*meas_cf )dΩ
    + ∫( q*(∇⋅u_int) )dΩ
    )

    if Dc == 2
      return v -> _liformX(v)
    elseif Dc == 3
      # in 3D, account for the boundary term from IBP
      Γ = BoundaryTriangulation(panel_model;tags=["bottom_boundary","top_boundary"])
      dΓ = Measure(Γ,degree)
      nΓ = get_normal_vector(Γ)
      boundary((v,q)) = ∫( (v⋅nΓ)*p_int )dΓ
      return v -> _liformX(v) - boundary(v)
    end
  end

  op = AffineFEOperator(biformX,get_liform(Dc),X,Y)
  uh,ph = solve(ls,op)

  uh_proj = covarient_basis_cf ⋅ (1/meas_cf*uh)

  _e = u_cf - uh
  e_u =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*(1/meas_cf) )dΩ_error))

  _e = h_cf - ph
  e_p = sqrt(sum(∫( (_e*_e)*meas_cf )dΩ_error))

  if return_vtk
    panel_cfs = [h_cf, u_proj_cf, ph, uh_proj, h_cf-ph, u_proj_cf-uh_proj ]
    labels = ["p","u_proj", "ph", "uh_proj", "ep","eu"]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$(p_fe)_D$Dc",
            cellfields=cellfields,append=false,geo_map=geo_map_func(Ω_panel))
  end

  return e_u,e_p,false
end

################################################################################
#### Launch wave equation -- on gadi
################################################################################
function launch_wave_equation(ranks,Dc,n_ref,p_fe::Int,dir::String,return_vtk=1)

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("Wave equation: Dc = $Dc")

  dir_convergence = dir*"/convergence"
  (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

  # ensure no MPI task tries to generate the file before the main MPI task has
  # created the folder
  PartitionedArrays.barrier(ranks)

  ζ = 0.0
  h = panel_to_cartesian(h₀(ζ))
  vX = panel_to_cartesian(tangent_vec(u₀(ζ)))

  omodel = if Dc == 2
    ParametricOctreeDistributedDiscreteModel(ranks;
    num_initial_uniform_refinements=n_ref)
  elseif Dc == 3
    Parametric3DOctreeDistributedDiscreteModel(ranks;
        num_horizontal_uniform_refinements=n_ref,
        num_vertical_uniform_refinements=n_ref);
  end

  panel_model = omodel.parametric_dmodel


  GridapPETSc.Init()
  ls = PETScLinearSolver(petsc_mumps_setup)

  e_u,e_p, = wave_solver(panel_model,p_fe,dir,h,vX,ls,Bool(return_vtk))

  i_am_main(ranks) && println("eu = $e_u, e_p = $e_p")

  ## convergence output for DrWatson
  n = nc(panel_model)
  dxx = dx(panel_model)
  output = @strdict e_u e_p n dxx p_fe n_ref Dc
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("wave_equation_nref$(n_ref)_p$(p_fe)_D$Dc.jld2")), output)


  GridapPETSc.Finalize()
  GridapPETSc.gridap_petsc_gc()

  i_am_main(ranks) && println("--DONE--")

end



################################################################################
#### Auto convergence test
################################################################################
function main(models::AbstractArray)
  h = panel_to_cartesian(h₀(0.0))
  vX = panel_to_cartesian(tangent_vec(u₀(0.0)))

  ls = LUSolver()
  dir = @__DIR__
  ps = [1,2]
  p_convergence_auto_test(ps,models,wave_solver,dir,h,vX,ls)
end

function main(distribute,nprocs;)
  ranks = distribute(LinearIndices((nprocs,)))

  n_ref_lvls = 4

  ## Distributed model: 2D
  models = get_distributed_refined_models(ranks,nprocs,n_ref_lvls)
  main(models)

  ### P4test model: 2D
  models = get_octree_refined_models(ranks,n_ref_lvls)
  main(models)

  ### P4test model: 3D
  models = get_3D_octree_refined_models(ranks,n_ref_lvls-1)
  main(models)

end
