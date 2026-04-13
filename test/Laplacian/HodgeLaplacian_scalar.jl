"""
solve scalar laplacian in mixed form
u + ∇ᵧ(φ) = 0
∇ᵧ⋅u = f
where f = -Δφ
"""

module HodgeLaplacianScalarTests

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
  function _f(α)
    xyz = ForwardMap(p)(α)
    θϕr   = xyz2θϕr(xyz)
    sin(θϕr[2])
  end
end

function hodge_laplacian_scalar(panel_model,
  p_fe::Int,dir::String,f::Function,ls=LUSolver(),return_vtk=false;
  _i_am_main=true)

  Dc = num_cell_dims(panel_model)
  lvl = nref(panel_model)
 _i_am_main && println("p_fe = $(p_fe); nref = $lvl; Dc = $Dc")

  degree = 4*(p_fe+1)
  if p_fe == 0
    degree = 10
  end
  @check degree > 0 "Zero quad!!"

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)
  dΩ_error = Measure(Ω_panel,2*degree)

  # FE spaces
  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  if Dc == 2
   _i_am_main && println("zeromean constraint in 2D ")
    Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2, constraint=:zeromean)
  end
  P = TrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  # metric information
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  # manufactured RHS
  f_panel_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  sigma_cf = panelwise_cellfield(sgrad(f),Ω_panel,panel_ids)
  slap_panel_cf =  panelwise_cellfield(surflap(f),Ω_panel,panel_ids)
  rhs = -slap_panel_cf
  f_int = interpolate(f_panel_cf,P)

  biform_u((u,p),(v,q)) = ∫( (u⋅ (metric_cf⋅v))*(1/meas_cf) )dΩ - ∫( p*(∇⋅v) )dΩ
  biform_p((u,p),(v,q)) = ∫( q*(∇⋅u) )dΩ

  biformX((u,p),(v,q)) = biform_u((u,p),(v,q)) + biform_p((u,p),(v,q))


  # manufacture rhs functions
  function get_liform(Dc::Int)

    # the manufactured solution
    _liformX((v,q)) = ∫( (rhs*q)*meas_cf )dΩ

    if Dc == 2
      return v -> _liformX(v)
    elseif Dc == 3
      # in 3D, account for the boundary term from IBP
      Γ = BoundaryTriangulation(panel_model;tags=["bottom_boundary","top_boundary"])
      dΓ = Measure(Γ,degree)
      nΓ = get_normal_vector(Γ)
      boundary((v,q)) = ∫( -f_int*(v⋅nΓ) )dΓ
      return v -> _liformX(v) + boundary(v)
    end
  end


  op = AffineFEOperator(biformX,get_liform(Dc),X,Y)
  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  xh = FEFunction(X,x)
  uh,ph = xh

  _e = f_panel_cf - ph
  el2_p = sqrt(sum(∫( (_e*_e)*meas_cf  )dΩ_error))

  _e = (covarient_basis_cf⋅(1/meas_cf*uh)) - (- sigma_cf ) ### u = -∇p
  el2_u = sqrt(sum(∫( (_e⋅_e)*meas_cf  )dΩ_error))

 _i_am_main && println("eu = $(el2_u), es = $(el2_p)")

  if return_vtk
    cellfields =  ["u"=> -sigma_cf ,
    "uh"=>covarient_basis_cf⋅(1/meas_cf*uh),
    "eu"=> (covarient_basis_cf⋅(1/meas_cf*uh)) - (-sigma_cf),
    "ph"=>ph, "p"=>f_panel_cf, "e"=>ph-f_panel_cf
                  ]
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",
            cellfields=cellfields,
            append=false,geo_map= geo_map_func(Ω_panel))
  end


  return el2_u, el2_p, false

end


################################################################################
#### Launch wave equation -- on gadi
################################################################################
# function launch_hodge_laplacian(ranks,n_ref,p_fe::Int,dir::String,return_vtk=1)

#   i_am_main(ranks) && println("--START--")
#   i_am_main(ranks) && println("Hodge Laplacian: scalar")

#   (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

#   dir_convergence = dir*"/convergence"
#   (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

#   # ensure no MPI task tries to generate the file before the main MPI task has
#   # created the folder
#   PartitionedArrays.barrier(ranks)

#   octree3_model = Parametric3DOctreeDistributedDiscreteModel(ranks;
#     num_horizontal_uniform_refinements=n_ref,
#     num_vertical_uniform_refinements=n_ref);
#   panel_model = octree3_model.parametric_dmodel

#   GridapPETSc.Init()
#   ls = PETScLinearSolver(petsc_mumps_setup)
#   # ls = LUSolver()

#   e_u, e_p, = hodge_laplacian_scalar(panel_model,p_fe,dir,fX,ls,Bool(return_vtk);_i_am_main=i_am_main(ranks))

#   ### convergence output for DrWatson
#   n = nc(panel_model)
#   dxx = dx(panel_model)
#   output = @strdict e_u e_p n dxx p_fe n_ref
#   i_am_main(ranks) && safesave(datadir(dir_convergence, ("hodge_laplacian_scalar_nref$(n_ref)_p$p_fe.jld2")), output)


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
  p_convergence_auto_test(ps,models,hodge_laplacian_scalar,dir,fX,ls;_i_am_main=_i_am_main)
end

function main(distribute,nprocs;)
  ranks = distribute(LinearIndices((nprocs,)))

  n_ref_lvls = 4

  ## Distributed model: 2D
  models = get_distributed_refined_models(ranks,nprocs,n_ref_lvls)
  main(models;_i_am_main=i_am_main(ranks))

  ### P4test model: 2D
  models = get_octree_refined_models(ranks,n_ref_lvls)
  main(models;_i_am_main=i_am_main(ranks))

  ### P4test model: 3D
  models = get_3D_octree_refined_models(ranks,n_ref_lvls-1)
  main(models;ps=[1],_i_am_main=i_am_main(ranks))

end



end # module
