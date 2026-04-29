""" Poisson problem using Laplace-Beltrami operator
u + Δᵧ(u) = f
Need to remove the kernal via zeromean FE space
"""

module LaplaceBeltramiTests

using Gridap
using Gridap.Helpers
using Gridap.Algebra
using GridapGeosciences
using GridapP4est
using Test

# using DrWatson
# using MPI
# using PartitionedArrays

# using GridapPETSc
# function petsc_mumps_setup(ksp)
#   pc       = Ref{GridapPETSc.PETSC.PC}()
#   mumpsmat = Ref{GridapPETSc.PETSC.Mat}()
#   @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPPREONLY)
#   @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
#   @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCLU)
#   @check_error_code GridapPETSc.PETSC.PCFactorSetMatSolverType(pc[],GridapPETSc.PETSC.MATSOLVERMUMPS)
#   @check_error_code GridapPETSc.PETSC.PCFactorSetUpMatSolverType(pc[])
#   @check_error_code GridapPETSc.PETSC.PCFactorGetMatrix(pc[],mumpsmat)
#   @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  4, 1)
#   @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 28, 2)
#   @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 29, 1)
#   # @check_error_code GridapPETSc.PETSC.MatMumpsSetCntl(mumpsmat[],  1, 0.00001)
#   @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
# end


function fX(p)
  function _f(αβ)
    x = ForwardMap(p)(αβ)
    x[1]*x[2]*x[3]
  end
end


function laplace_beltrami_solver(panel_model,
  p_fe::Int,dir::String,f::Function,ls=LUSolver(),return_vtk=false;
  _i_am_main=true)


  Dc = num_cell_dims(panel_model)
  lvl = nref(panel_model)

  _i_am_main && println("nref = $lvl; p_fe = $p_fe; Dc = $Dc")

  degree = 6*(p_fe+1)
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)
  dΩ_error = Measure(Ω_panel,2*degree)

  V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1, constraint=:zeromean)
  U = TrialFESpace(V)

  f_panel_cf = ParametricCellField(f,Ω_panel,panel_ids)
  inv_metric_cf = ParametricCellField(inv_metric,Ω_panel,panel_ids)
  meas_cf = ParametricCellField(sqrtg,Ω_panel,panel_ids)
  slap_panel_cf =  ParametricCellField(surflap(f),Ω_panel,panel_ids)

  @check sum(∫(f_panel_cf*meas_cf)dΩ) < 1e-14 "Function must be zero mean to solve with zeromean FE space!"

  rhs_cf = - slap_panel_cf

  poisson_biform(u,v) =  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_cf )dΩ

  # manufacture rhs functions
  function get_liform(Dc::Int)

    # the manufactured solution
    poisson_liform(v) = ∫(  (rhs_cf*v)*meas_cf )dΩ

    if Dc == 2
      return v -> poisson_liform(v)
    elseif Dc == 3
      # in 3D, account for the boundary term from IBP
      Γ = BoundaryTriangulation(panel_model;tags=["bottom_boundary","top_boundary"])
      dΓ = Measure(Γ,degree)
      nΓ = get_normal_vector(Γ)
      f_int = interpolate(f_panel_cf,U)
      boundary(v) = ∫( ( (inv_metric_cf⋅gradient(f_int) )⋅nΓ)*v*meas_cf )dΓ
      return v -> poisson_liform(v) + boundary(v)
    end
  end

  op = AffineFEOperator(poisson_biform,get_liform(Dc),U,V)

  # uh = solve(ls,op)

  ## for pvectors, the ghost may not be in the prange of the get_matrix
  ## This causes issues with GridapSolvers Krylov solvers, in the allocation of x
  ## To avoid, allocate x based on the domain of A
  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  uh = FEFunction(U,x)

  _e = f_panel_cf - uh
  e = sqrt(sum(∫( (_e*_e)*meas_cf )dΩ_error))

  if return_vtk
    panel_cfs = [f_panel_cf,uh,f_panel_cf-uh]
    labels = ["u","uh","eu"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk_with_cell_geomap(geo_map_func(Ω_panel),Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",
        cellfields=cellfields,append=false)
  end

  ### convergence output for DrWatson


  return e, false,false
end



################################################################################
#### Launch launch_laplace_beltrami -- on gadi
################################################################################
# function launch_laplace_beltrami(ranks,Dc,n_ref,p_fe::Int,dir::String,return_vtk=1)

#   i_am_main(ranks) && println("--START--")
#   i_am_main(ranks) && println("Laplace Beltrami: Dc = $Dc")

#   (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

#   dir_convergence = dir*"/convergence"
#   (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

#   # ensure no MPI task tries to generate the file before the main MPI task has
#   # created the folder
#   PartitionedArrays.barrier(ranks)

#   omodel = if Dc == 2
#     ParametricOctreeDistributedDiscreteModel(ranks, radius;
#     num_initial_uniform_refinements=n_ref)
#   elseif Dc == 3
#     Parametric3DOctreeDistributedDiscreteModel(ranks,radius,thickness;
#         num_horizontal_uniform_refinements=n_ref,
#         num_vertical_uniform_refinements=n_ref);
#   end

#   panel_model = omodel.parametric_dmodel


#   GridapPETSc.Init()
#   ls = PETScLinearSolver(petsc_mumps_setup)
#   f = fX
#   e_u,  = laplace_beltrami_solver(panel_model,p_fe,dir,f,ls,Bool(return_vtk);_i_am_main=i_am_main(ranks))

#   i_am_main(ranks) && println("eu = $e_u")

#   ## convergence output for DrWatson
#   n = nc(panel_model)
#   dxx = dx(panel_model)
#   output = @strdict e_u n dxx p_fe n_ref Dc
#   i_am_main(ranks) && safesave(datadir(dir_convergence, ("laplace_beltrami_nref$(n_ref)_p$(p_fe)_D$Dc.jld2")), output)

#   GridapPETSc.Finalize()
#   GridapPETSc.gridap_petsc_gc()

#   i_am_main(ranks) && println("--DONE--")

# end



################################################################################
#### Auto convergence test
################################################################################
function main(models::AbstractArray;ps=[2],_i_am_main=true)

  ls = LUSolver()
  dir = @__DIR__
  p_convergence_auto_test(ps,models,laplace_beltrami_solver,dir,fX,ls;_i_am_main=_i_am_main)
end

# function main(distribute,nprocs;)
#   ranks = distribute(LinearIndices((nprocs,)))

#   n_ref_lvls = 4
#   radius = 1
# thickness = 0.19
#   ## Distributed model: 2D
#   models = get_distributed_refined_models(ranks,nprocs,n_ref_lvls,radius)
#   main(models;_i_am_main=i_am_main(ranks))

#   ### P4test model: 2D
#   models = get_octree_refined_models(ranks,n_ref_lvls,radius)
#   main(models;_i_am_main=i_am_main(ranks))

#   ### P4test model: 3D
#   models = get_3D_octree_refined_models(ranks,n_ref_lvls-1,radius,thickness)
#   main(models;_i_am_main=i_am_main(ranks))

# end



end # module
