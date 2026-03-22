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

include("../convergence_tools.jl")

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

  octree3_model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
    num_horizontal_uniform_refinements=n_ref,
    num_vertical_uniform_refinements=n_ref);
  panel_model = octree3_model.parametric_dmodel

  GridapPETSc.Init()
  ls = PETScLinearSolver(petsc_mumps_setup)

  hodge_laplacian(panel_model,p_fe,dir,uX,ls,Bool(return_vtk))

  GridapPETSc.Finalize()
  GridapPETSc.gridap_petsc_gc()

  i_am_main(ranks) && println("--DONE--")

end

function hodge_laplacian(
  panel_model::GridapDistributed.GenericDistributedDiscreteModel{3,3},
  p_fe::Int,dir::String,uX::Function,ls=LUSolver(),return_vtk=false)

  ranks = get_ranks(panel_model)
  lvl_h = nref(nc_horizontal(panel_model))
  lvl_v = nref(nc_vertical(panel_model))
  i_am_main(ranks) && println("nref_h = $lvl_h; nref_v = $lvl_v; p_fe = $p_fe")


  degree = 5*(p_fe + 1)
  if p_fe == 0
    degree = 10
  end
  @check degree > 0 "Zero quad!!"

  ## finite element solver
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*degree)
  Ω_error = Triangulation(panel_model)
  dΩ_error = Measure(Ω_error,4*degree)

  tags = ["top_boundary", "bottom_boundary"]
  Γ = BoundaryTriangulation(panel_model,tags=tags)
  dΓ = Measure(Γ,2*degree)
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


  cellfields = ["curlu"=>ccurlu_cov_cf,
                "u"=>covarient_basis_cf ⋅ (inv_metric_cf⋅u_cov_cf),
                "un"=>un_cf,
                "curlu_cross"=>covarient_basis_cf ⋅ (inv_metric_cf⋅curlu_cross),
                "sigma"=>-sdiv_cf,
                "rhs"=>rhs_cov_cf
                ]
  writevtk(Ω_panel,dir*"/sol",
          cellfields=cellfields,
          append=false,geo_map= geo_map_func(Ω_panel))


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
  sh,uh = solve(ls,op)


  # _e = sigma_cf - sh
  _e = sigma_int - sh
  el2_s = sqrt(sum(∫( (_e*_e)*meas_cf  )dΩ_error))

  # _e = (inv_metric_cf⋅uh) - (inv_metric_cf⋅u_cov_cf)
  _e = (inv_metric_cf⋅uh) - (inv_metric_cf⋅u_int)
  el2_u = sqrt(sum(∫( (_e⋅(metric_cf ⋅_e))*meas_cf  )dΩ_error))

  if return_vtk
    cellfields =  ["u"=>covarient_basis_cf ⋅ (inv_metric_cf⋅u_cov_cf),
    "uh"=>covarient_basis_cf ⋅ (inv_metric_cf⋅uh),
    "eu"=>covarient_basis_cf ⋅ (inv_metric_cf⋅uh)-covarient_basis_cf ⋅ (inv_metric_cf⋅u_cov_cf),
    "sh"=>sh, "s"=>sigma_cf, "e"=>sh-sigma_cf
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
  output = @strdict el2_u el2_s n n_h n_v dxx dxH dxV p_fe lvl_h lvl_v
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("hodge_laplacian_nrefh$(lvl_h)_nrefv$(lvl_v)_p$p_fe.jld2")), output)

  return el2_u, el2_s, false

end


# MPI.Init()
# ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

# dir = datadir("HodgeLaplacianConvergence")
# (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

# n_ref_lvls = 3
# ps = [0,1]
# models  = get_3D_octree_refined_models(ranks,n_ref_lvls)
# ls = LUSolver()

# p_convergence_test(ranks,ps,models,hodge_laplacian,dir,uX,ls,true)
# i_am_main(ranks) && plot_convergence_from_saved(dir,"convergence",["u:", "s:"])
