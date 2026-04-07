""" Helmholtz problem
u + Δᵧ(u) = f
"""

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using GridapGeosciences
using Test

using MPI
using PartitionedArrays

function fX(p)
  function _f(αβ)
    x = ForwardMap(p)(αβ)
    x[1]*x[2]*x[3]
  end
end

function helmholtz_solver(panel_model,
  p_fe::Int,dir::String,f::Function,ls=LUSolver(),return_vtk=false)

  ranks = get_ranks(panel_model)
  Dc = num_cell_dims(panel_model)
  lvl = nref(panel_model)

  i_am_main(ranks) && println("nref = $lvl; p_fe = $p_fe; Dc = $Dc")

  degree = 6*(p_fe+1)
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)
  dΩ_error = Measure(Ω_panel,2*degree)

  V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  U = TrialFESpace(V)

  f_panel_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  slap_panel_cf = panelwise_cellfield(surflap(f),Ω_panel,panel_ids)

  rhs_cf = f_panel_cf + slap_panel_cf

  helmholtz_biform(u,v) = ∫(u*v*meas_cf)dΩ -  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_cf )dΩ

  # manufacture rhs functions
  function get_liform(Dc::Int)

    # the manufactured solution
    helmholtz_liform(v) = ∫( (rhs_cf*v)*meas_cf )dΩ

    if Dc == 2
      return v -> helmholtz_liform(v)
    elseif Dc == 3
      # in 3D, account for the boundary term from IBP
      Γ = BoundaryTriangulation(panel_model;tags=["bottom_boundary","top_boundary"])
      dΓ = Measure(Γ,degree)
      nΓ = get_normal_vector(Γ)
      f_int = interpolate(f_panel_cf,U)
      boundary(v) = ∫( ( (inv_metric_cf⋅gradient(f_int) )⋅nΓ)*v*meas_cf )dΓ
      return v -> helmholtz_liform(v) + boundary(v)
    end
  end

  op = AffineFEOperator(helmholtz_biform,get_liform(Dc),U,V)


  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = Gridap.Algebra.allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  uh = FEFunction(U,x)

  _e = f_panel_cf - uh
  e = sqrt(sum(∫( (_e*_e)*meas_cf )dΩ_error))

  if return_vtk
    panel_cfs = [f_panel_cf,uh,f_panel_cf-uh]
    labels = ["u","uh","eu"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",cellfields=cellfields,append=false,geo_map=geo_map_func(Ω_panel))
  end

  return e, false, false
end



################################################################################
#### Launch Helmholtz -- on gadi
################################################################################
function launch_helmholtz(ranks,Dc,n_ref,p_fe::Int,dir::String,return_vtk=1)

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("Helmholtz: Dc = $Dc")

  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  dir_convergence = dir*"/convergence"
  (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

  # ensure no MPI task tries to generate the file before the main MPI task has
  # created the folder
  PartitionedArrays.barrier(ranks)

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
  f = fX
  e_u,  = helmholtz_solver(panel_model,p_fe,dir,f,ls,Bool(return_vtk))

  i_am_main(ranks) && println("eu = $e_u")

  ## convergence output for DrWatson
  n = nc(panel_model)
  dxx = dx(panel_model)
  output = @strdict e_u n dxx p_fe n_ref Dc
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("helmholtz_nref$(n_ref)_p$(p_fe)_D$Dc.jld2")), output)

  GridapPETSc.Finalize()
  GridapPETSc.gridap_petsc_gc()

  i_am_main(ranks) && println("--DONE--")

end



################################################################################
#### Auto convergence test
################################################################################
function main(models::AbstractArray)
  f = fX

  ls = LUSolver()
  dir = @__DIR__
  ps = [1,2,3]

  p_convergence_auto_test(ps,models,helmholtz_solver,dir,f,ls)
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
