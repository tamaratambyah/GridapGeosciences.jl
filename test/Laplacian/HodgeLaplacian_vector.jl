"""
solve vector laplacian in mixed form
σ + ∇ᵧ⋅ϕ  = 0
∇ᵧ × (∇ᵧ × ϕ) + ∇ᵧσ = f
where f = -Δϕ
"""

module HodgeLaplacianVectorTests

using MPI
using PartitionedArrays

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra
using Gridap.ReferenceFEs, Gridap.Polynomials, Gridap.CellData

using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Test
using GridapPETSc
using LinearAlgebra

inv_jacobian(p) = x -> inv(forward_jacobian_3D(p)(x))
contra_v_3D(vecX::Function,p::Int) = x -> inv_jacobian(p)(x) ⋅ vecX(p)(x)
contra_v_3D(vecX::Function) = p -> contra_v_3D(vecX,p)

transpose_jacobian(p) = x -> transpose(forward_jacobian_3D(p)(x))
inv_tranpose_jacobian(p) = x -> inv(transpose_jacobian(p)(x))

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

function uX(p)
  function _u(γαβ)
    xyz = forward_map_3D(p)(γαβ)
    # VectorValue(-xyz[2],xyz[1],0.0)

    r = sqrt(xyz[1]^2 + xyz[2]^2 + xyz[3]^2)
    f = 2.0*xyz[3]/r
    n = normal_vec(xyz)
    f*n
  end
end


function launch_hodge_laplacian(ranks,n_ref,p_fe::Int,dir::String,return_vtk=1)

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("Hodge")

  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  ### convergence output for DrWatson
  dir_convergence = dir*"/convergence"
  (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

  octree3_model = Parametric3DOctreeDistributedDiscreteModel(ranks;
    num_horizontal_uniform_refinements=n_ref,
    num_vertical_uniform_refinements=n_ref);
  panel_model = octree3_model.parametric_dmodel

  GridapPETSc.Init()
  ls = PETScLinearSolver(petsc_mumps_setup)

  e_u, e_s, =  hodge_laplacian_vector(panel_model,p_fe,dir,uX,ls,Bool(return_vtk);_i_am_main=i_am_main(ranks))

  ### convergence output for DrWatson
  n = nc(panel_model)
  dxx = dx(panel_model)
  output = @strdict e_u e_s n dxx p_fe n_ref
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("hodge_laplacian_vector_nref$(n_ref)_p$p_fe.jld2")), output)


  GridapPETSc.Finalize()
  GridapPETSc.gridap_petsc_gc()

  i_am_main(ranks) && println("--DONE--")

end

function hodge_laplacian_vector(
  panel_model::GridapDistributed.GenericDistributedDiscreteModel{3,3},
  p_fe::Int,dir::String,uX::Function,ls=LUSolver(),return_vtk=false;
  _i_am_main=true)

  Dc = num_cell_dims(panel_model)
  lvl = nref(panel_model)
 _i_am_main && println("p_fe = $(p_fe); nref = $lvl; Dc = $Dc")

  # degree = 30
  degree = 5*(p_fe + 1)
  if p_fe == 0
    degree = 10
  end
  @check degree > 0 "Zero quad!!"

  ## finite element solver
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)
  Ω_error = Triangulation(panel_model)
  dΩ_error = Measure(Ω_error,2*degree)

  tags = ["top_boundary", "bottom_boundary"]
  Γ = BoundaryTriangulation(panel_model,tags=tags)
  dΓ = Measure(Γ,degree)
  nΓ = get_normal_vector(Γ)

  ## metric information
  inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)



  # covarient components of u
  function ucov(p)
    function _u(γαβ)
      u = uX(p)(γαβ)
      J = transpose_jacobian(p)(γαβ)
      J⋅u
    end
  end

  # u ⋅ n on the surface
  function unX(p)
    function _u(γαβ)
      xyz = forward_map_3D(p)(γαβ)
      u = uX(p)(γαβ)
      n = normal_vec(xyz)
      u⋅n
    end
  end

  ### Curl of covariant components of u
  curlu(p,x) = curl(ucov(p))(x)
  curlu(p) = x -> curlu(p,x)
  _curlu(p) = curl(ucov(p))

  ### Covariant components of surfcurl u
  wcov(p,x) = 1/sqrtg(p,x)*metric(p,x)⋅curlu(p,x)
  wcov(p) = x -> wcov(p,x)
  curlw(p,x) = curl(wcov(p))(x)
  curlw(p) = x -> curlw(p,x)

  #### Covariant component of surfcurl surfcurl u
  curlw_cov(p) = x -> 1/sqrtg(p,x)*metric(p,x)⋅curlw(p,x)

  # area measure
  _area_meas(p) = x->  forward_jacobian_3D(p,x) ⋅ (inv_metric(p,x) ⋅ VectorValue(1,0,0))
  area_meas(p) = x-> norm(_area_meas(p)(x))

  #### Covariant componetsn of (surfcurl u)× surfnormal
  wcrossk_cov(p) = x -> 1/sqrtg(p,x) * metric(p,x)⋅(wcov(p,x) × (VectorValue(1,0,0)/area_meas(p)(x)) )


  # _t(p) = x -> sqrtg(p)(x)* contra_v_3D(uX,p)(x)
  # _k(p) = x -> 1/sqrtg(p)(x) * ( divergence(_t(p) )(x) )
  # _l(p) = gradient(_k(p))
  # rhs(p) = x-> curlw_cov(p)(x) - _l(p)(x)

  ## surface divergence
  _sdiv_u(p) = x -> sqrtg(p)(x)* contra_v_3D(uX,p)(x)
  sdiv_u(p) = x -> 1/sqrtg(p)(x) * ( divergence(_sdiv_u(p) )(x) )

  ### covariant components of surfgrad(surfdiv u)
  graddiv_cov(p) = gradient(sdiv_u(p))

  ### Rhs function
  rhs(p) = x-> curlw_cov(p)(x) - graddiv_cov(p)(x)
  rhs_cov_cf = panelwise_cellfield(rhs,Ω_panel,panel_ids)

  u_cov_cf = panelwise_cellfield(ucov,Ω_panel,panel_ids)
  ccurlu_cov_cf = panelwise_cellfield(curlw_cov,Ω_panel,panel_ids)
  un_cf = panelwise_cellfield(unX,Ω_panel,panel_ids)
  curlu_cross = panelwise_cellfield(wcrossk_cov,Ω_panel,panel_ids)
  curlu_cf = panelwise_cellfield(curlu,Ω_panel,panel_ids)

  sdiv_cf =  panelwise_cellfield(surfdiv(contra_v_3D(uX)),Ω_panel,panel_ids)
  sigma_cf = -sdiv_cf


  # cellfields = ["curlu"=>ccurlu_cov_cf,
  #               "u"=>covarient_basis_cf ⋅ (inv_metric_cf⋅u_cov_cf),
  #               "un"=>un_cf,
  #               "curlu_cross"=>covarient_basis_cf ⋅ (inv_metric_cf⋅curlu_cross),
  #               "sigma"=>-sdiv_cf,
  #               "rhs"=>rhs_cov_cf
  #               ]
  # writevtk(Ω_panel,dir*"/sol",
  #         cellfields=cellfields,
  #         append=false,geo_map= geo_map_func(Ω_panel))


  ## FE spaces
  T = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1)
  S = TrialFESpace(T)

  R = TestFESpace(Ω_panel, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl)
  H = TrialFESpace(R)

  sigma_int = interpolate(sigma_cf,S)
  u_int = interpolate(u_cov_cf,H)

  ### Multifield
  X = MultiFieldFESpace([S,H])
  Y = MultiFieldFESpace([T,R])

  biform_x((s,u),(t,v)) = (
                  ∫( (s*t)*meas_cf  )dΩ
                - ∫( ∇(t)⋅(inv_metric_cf⋅u)*meas_cf  )dΩ
                + ∫( curl(u)⋅(metric_cf⋅curl(v))*(1/meas_cf) )dΩ
                + ∫( gradient(s)⋅(inv_metric_cf⋅v)*meas_cf )dΩ
                  )
  liform_x((t,v)) = (
                ∫( rhs_cov_cf⋅(inv_metric_cf⋅v)*meas_cf  )dΩ
                + ∫( v⋅( ( metric_cf⋅curlu_cf )×nΓ    )*(1/meas_cf)     )dΓ
                - ∫(( t*(u_cov_cf⋅(inv_metric_cf⋅nΓ)) )*(meas_cf)  )dΓ
                  )


  op = AffineFEOperator(biform_x,liform_x,X,Y)
  xh = solve(ls,op)
  sh,uh = xh

  # final_dir = dir*"/final_solution"
  # # ensure no MPI task tries to generate the file before the main MPI task has
  # # created the folder
  # PartitionedArrays.barrier(ranks)
  # psave(final_dir*"/sol",xh)

  # _e = sigma_cf - sh
  _e = sigma_int - sh
  el2_s = sqrt(sum(∫( (_e*_e)*meas_cf  )dΩ_error))

  # _e = (inv_metric_cf⋅uh) - (inv_metric_cf⋅u_cov_cf)
  _e = (inv_metric_cf⋅uh) - (inv_metric_cf⋅u_int)
  el2_u = sqrt(sum(∫( (_e⋅(metric_cf ⋅_e))*meas_cf  )dΩ_error))

 _i_am_main && println("eu = $(el2_u), es = $(el2_s)")

  if return_vtk
    cellfields =  ["u"=>covarient_basis_cf ⋅ (inv_metric_cf⋅u_cov_cf),
    "uh"=>covarient_basis_cf ⋅ (inv_metric_cf⋅uh),
    "eu"=>covarient_basis_cf ⋅ (inv_metric_cf⋅uh)-covarient_basis_cf ⋅ (inv_metric_cf⋅u_cov_cf),
    "sh"=>sh, "s"=>sigma_cf, "e"=>sh-sigma_cf
                  ]
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",
            cellfields=cellfields,
            append=false,geo_map= geo_map_func(Ω_panel))
  end


  return el2_u, el2_s, false

end


################################################################################
#### Auto convergence test -- only 3D models for p=[0,1]
################################################################################
function main(models::AbstractArray;ps = [1],_i_am_main=true)

  ls = LUSolver()
  dir = @__DIR__
  p_convergence_auto_test(ps,models,hodge_laplacian_vector,dir,uX,ls;_i_am_main=_i_am_main)
end

function main(distribute,nprocs;)
  ranks = distribute(LinearIndices((nprocs,)))

  n_ref_lvls = 3

  ### P4test model: 3D
  models = get_3D_octree_refined_models(ranks,n_ref_lvls)
  main(models;_i_am_main=i_am_main(ranks))

end



end # module
