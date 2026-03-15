using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed
using DrWatson
using Gridap.ReferenceFEs, Gridap.Polynomials

include("../../convergence_tools.jl")
# include("../../Geophysical/CurlConformingFESpacesFixes.jl")

## pullback 3D vector to 3D chart
inv_jacobian(p) = x -> inv(forward_jacobian_3D(p)(x))
contra_v_3D(vecX::Function,p::Int) = αβ -> inv_jacobian(p)(αβ) ⋅ vecX(p)(αβ)
contra_v_3D(vecX::Function) = p -> contra_v_3D(vecX,p)

transpose_jacobian(p) = x -> transpose(forward_jacobian_3D(p)(x))
inv_tranpose_jacobian(p) = x -> inv(transpose_jacobian(p)(x))
covar_v_3D(vecX::Function,p::Int) = αβ -> transpose_jacobian(p)(αβ) ⋅ vecX(p)(αβ)
covar_v_3D(vecX::Function) = p -> covar_v_3D(vecX,p)

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

function Gridap.ReferenceFEs._Nedelec_edge_values(p,et,order)
  println("my nedelec")

  # Reference facet
  dim1 = 1
  ep = Polytope{dim1}(p,1)

  # geomap from ref face to polytope faces
  egeomap = Gridap.ReferenceFEs._ref_face_to_faces_geomap(p,ep)

  # Compute integration points at all polynomial edges
  degree = (order)*2 + 2
  equad = Quadrature(ep,degree)
  cips = get_coordinates(equad)
  wips = get_weights(equad)


  c_eips, ecips, ewips = Gridap.ReferenceFEs._nfaces_evaluation_points_weights(p, egeomap, cips, wips)

  # Edge moments, i.e., M(Ei)_{ab} = q_RE^a(xgp_REi^b) w_Fi^b t_Ei ⋅ ()
  eshfs = MonomialBasis(et,ep,order)
  emoments = Gridap.ReferenceFEs._Nedelec_edge_moments(p, eshfs, c_eips, ecips, ewips)

  return ecips, emoments

end


function interpolation(panel_model::GridapDistributed.GenericDistributedDiscreteModel{3,3},
  p_fe::Int,dir::String,vecX::Function,conf,return_vtk)

  Dc = num_cell_dims(panel_model)
  @check Dc == 3

  lvl = nref(nc_horizontal(panel_model))
  println("p_fe = $(p_fe); nref = $lvl; Dc = $Dc")

  degree = 4*p_fe
  if p_fe == 0
    degree = 8
  end

  println("p_fe = $(p_fe); degree = $(degree)")

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)
  panel_ids = get_panel_ids(panel_model)

  inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  vec_cov_cf = panelwise_cellfield(covar_v_3D(vecX),Ω_panel,panel_ids)
  vec_proj_cf = covarient_basis_cf⋅ ( inv_metric_cf ⋅ vec_cov_cf)


  reffe = ReferenceFE(nedelec,Float64,p_fe)
  tags = ["top_boundary", "bottom_boundary"]
  V = TestFESpace(panel_model, reffe; conformity=:Hcurl,dirichlet_tags=tags)
  U = TrialFESpace(V,vec_cov_cf)

  ## interpolation
  vec_cov_h = interpolate(vec_cov_cf,U)

  _e = ( inv_metric_cf ⋅ vec_cov_h) - ( inv_metric_cf ⋅ vec_cov_cf)
  el2_interp =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ))


  ## L2 projection
  a(u,v) = ∫( u⋅( inv_metric_cf⋅v)*meas_cf )dΩ
  l(v) = ∫( vec_cov_cf⋅(inv_metric_cf⋅v)*meas_cf )dΩ
  op = AffineFEOperator(a,l,U,V)
  vec_cov_h = solve(LUSolver(),op)

  vec_proj_h = covarient_basis_cf ⋅ ( inv_metric_cf ⋅ vec_cov_h)

  _e = ( inv_metric_cf ⋅ vec_cov_h) - ( inv_metric_cf ⋅ vec_cov_cf)
  el2_proj =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ))

  if return_vtk
    cellfields=["uproj"=> vec_proj_cf, "u_cov"=>vec_cov_cf,
    "uprojh"=>vec_proj_h, "u_covh"=>vec_cov_h,
    "eproj"=>vec_proj_cf-vec_proj_h, "e"=>vec_cov_cf-vec_cov_h]
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$(p_fe)_"*String(conf),cellfields=cellfields,
          append=false,geo_map=latlon_geo_map_func(Ω_panel))
  end

  return el2_interp, el2_proj, false

end


function vecX(p)
  function f(γαβ)
    xyz = forward_map_3D(p)(γαβ)
    # VectorValue(xyz[1],xyz[2],xyz[3])
    # VectorValue(xyz[2], 0.0, 0.0)
    VectorValue(0.0,xyz[3],xyz[1]^2)
  end
end

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))


n_ref_lvls = 3
ps = [0,1]

dir = datadir("InterpolationConvergence")
!isdir(dir) && mkdir(dir)

models  = get_3D_octree_refined_models(ranks,n_ref_lvls)
conf = :Hcurl
_dir = dir*"/vector_func_"*String(conf)
!isdir(_dir) && mkdir(_dir)
p_convergence_test(ranks,ps,models,interpolation,_dir,vecX,conf,true)
plot_convergence_from_saved(_dir,"convergence",["Interp:", "L2 proj:"])
