using MPI
using PartitionedArrays

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Test
using GridapPETSc

dir = datadir("SW_3D")
!isdir(dir) && mkdir(dir)

include("../convergence_tools.jl")
include("Williamson2Test_3D_testcase.jl")
# include("CurlConformingFESpacesFixes.jl")

inv_jacobian(p) = x -> inv(forward_jacobian_3D(p)(x))
contra_v_3D(vecX::Function,p::Int) = x -> inv_jacobian(p)(x) ⋅ vecX(p)(x)
contra_v_3D(vecX::Function) = p -> contra_v_3D(vecX,p)

transpose_jacobian(p) = x -> transpose(forward_jacobian_3D(p)(x))
inv_tranpose_jacobian(p) = x -> inv(transpose_jacobian(p)(x))
contravariant_basis_3D(p) = x -> inv_tranpose_jacobian(p)(x)

covar_v_3D(vecX::Function,p::Int) = x -> transpose_jacobian(p)(x) ⋅ vecX(p)(x)
covar_v_3D(vecX::Function) = p -> covar_v_3D(vecX,p)

contra_surfcurl(vec::Function,p::Int) = x -> 1/sqrtg(p,x) * (curl(covar_v_3D(vec,p))(x) )
contra_surfcurl(vec::Function) = p -> contra_surfcurl(vec,p)

surfcurl(vec::Function,p::Int) = x -> forward_jacobian_3D(p,x) ⋅ contra_surfcurl(vec,p)(x)
surfcurl(vec::Function) = p -> surfcurl(vec,p)

cov_surfcurl(vec::Function,p::Int) = x -> metric(p,x)⋅contra_surfcurl(vec,p)(x)
cov_surfcurl(vec::Function) = p -> cov_surfcurl(vec,p)


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

n_ref_lvls = 4

models  = get_3D_octree_refined_models(ranks,n_ref_lvls)
ls = LUSolver()
ps = [0,1]
p_convergence_test(ranks,ps,models,single_layer_shallow_water_vorticity,dir,ls,true)
plot_convergence_from_saved(dir,"convergence",["Err:", "Top:", "Bottom"])


function single_layer_shallow_water_vorticity(
  _panel_model::GridapDistributed.DistributedDiscreteModel{3,3},
  p_fe::Int,dir::String,
  ls=LUSolver(),return_vtk=false
  )

  ranks = get_ranks(_panel_model)

  lvl_h = nref(nc_horizontal(_panel_model))
  lvl_v = nref(nc_vertical(_panel_model))
  i_am_main(ranks) && println("nref_h = $lvl_h; nref_v = $lvl_v; p_fe = $p_fe")

  ## make new model with only 1 cell in vertical
  o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
        num_horizontal_uniform_refinements=lvl_h,
        num_vertical_uniform_refinements=0)
  panel_model = o3model.parametric_dmodel

  degree = 4*(p_fe + 1)

  if p_fe == 0
    degree = 8
  end

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

  ## cellfields
  u_contra_cf = panelwise_cellfield(contra_v_3D(u_vec_3D),Ω_panel,panel_ids)
  u_cov_cf = panelwise_cellfield(covar_v_3D(u_vec_3D),Ω_panel,panel_ids)
  f_cov_cf = panelwise_cellfield(covar_v_3D(f_vec_3D),Ω_panel,panel_ids)
  η_cov_cf = panelwise_cellfield(covar_v_3D(η_vec_3D),Ω_panel,panel_ids)
  n_cov = CellField(x->VectorValue(1,0,0),Ω_panel) ## same as nΓ
  h_cf = panelwise_cellfield(h_3D,Ω_panel,panel_ids)

  ## metric information
  inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  contravariant_basis_cf = panelwise_cellfield(contravariant_basis_3D,Ω_panel,panel_ids)
  jac_cf = panelwise_cellfield(forward_jacobian,Ω_panel,panel_ids)
  area_meas_cf = Operation(norm)(jac_cf⋅(inv_metric_cf ⋅nΓ) )

  ## manufacture rhs
  cov_scurl = panelwise_cellfield(cov_surfcurl(u_vec_3D),Ω_panel,panel_ids)
  rhs = η_cov_cf - cov_scurl - f_cov_cf


  ## FE spaces
  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv,dirichlet_tags=tags)
  U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))

  R = TestFESpace(panel_model, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl,dirichlet_tags=tags)
  H = TrialFESpace(R,VectorValue(0.0,0.0,0.0))

  # initial conditions
  u_contra_h = interpolate(u_contra_cf,U)
  η_cov_h = interpolate(η_cov_cf,H)

  function my_func(w,u)
    w[2]*u[3] - w[3]*u[2]
  end


  biformq(q,w) = ∫( h_cf*(q⋅(inv_metric_cf⋅w))*meas_cf )dΩ
  liformq(w) = (
                ∫( u_contra_h⋅( metric_cf⋅ curl(w) )  )dΩ
              - ∫( (( w × (metric_cf⋅ u_contra_h) )⋅nΓ)*area_meas_cf   )dΓ
                # - ∫( (my_func∘(w,metric_cf⋅ u_contra_h)) *area_meas_cf   )dΓ
              #  - ∫( (( w × u_cov_cf)⋅nΓ)*area_meas_cf   )dΓ
              #  - ∫( (my_func∘(w,u_cov_cf)) *area_meas_cf   )dΓ
              + ∫( (f_cov_cf⋅(inv_metric_cf ⋅ w))*meas_cf )dΩ
              + ∫(rhs⋅(inv_metric_cf⋅w)*meas_cf )dΩ
                )
  op = AffineFEOperator(biformq,liformq,H,R)
  qh = solve(ls,op)
  ηh = qh*h_cf

  ηh_ambient = covarient_basis_cf ⋅ (inv_metric_cf ⋅ ηh )
  η_ambient = covarient_basis_cf ⋅ (inv_metric_cf ⋅ η_cov_h )

  _e = ( inv_metric_cf ⋅ ηh) - ( inv_metric_cf ⋅  η_cov_h)
  el2_proj =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ_error))

  ### compute error on bottom boundary (inner shell )
  Γ_bottom = BoundaryTriangulation(panel_model,tags=["bottom_boundary"])
  dΓ_bottom = Measure(Γ_bottom,6)
  el2_bottom =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΓ_bottom))

  ### compute error on bottom boundary (outer shell )
  Γ_top = BoundaryTriangulation(panel_model,tags=["top_boundary"])
  dΓ_top = Measure(Γ_top,6)
  el2_top =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΓ_top))


  if return_vtk

    cellfields = ["ηambient"=>η_ambient,
                  "ηambient_h"=>ηh_ambient,
                  "eη_ambient"=>η_ambient-ηh_ambient,
                  "η_h"=>ηh ,
                  "eη"=>η_cov_h-ηh,  ]

    ### plot in 3D
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl_h)_p$p_fe",
            cellfields=cellfields,append=false,geo_map= geo_map_func(Ω_panel))

    ### Plot on bottom/top boundary
    writevtk(Γ_bottom,dir*"/ambient_model_bottom_nref$(lvl_h)_p$p_fe",
              cellfields=cellfields,append=false,geo_map=geo_map_func(get_panel_ids(Γ_bottom)))
    writevtk(Γ_top,dir*"/ambient_model_top_nref$(lvl_h)_p$p_fe",
            cellfields=cellfields,append=false,geo_map=geo_map_func(get_panel_ids(Γ_top)))
  end


  return el2_proj, el2_bottom, el2_top


end
