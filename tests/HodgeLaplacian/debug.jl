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

include("../convergence_tools.jl")

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


## increase quadrature in interpolation
function Gridap.ReferenceFEs._Nedelec_edge_values(p,et,order)
  println("Edges: increased quad")

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

function Gridap.ReferenceFEs._Nedelec_face_values(p,et,order)
  println("Faces: increased quad")

  # Reference facet
  @assert is_n_cube(p) "We are assuming that all n-faces of the same n-dim are the same."
  fp = Polytope{num_dims(p)-1}(p,1)

  # geomap from ref face to polytope faces
  fgeomap = Gridap.ReferenceFEs._ref_face_to_faces_geomap(p,fp)

  # Compute integration points at all polynomial edges
  degree = (order)*2 + 2
  fquad = Quadrature(fp,degree)
  fips = get_coordinates(fquad)
  wips = get_weights(fquad)

  c_fips, fcips, fwips = Gridap.ReferenceFEs._nfaces_evaluation_points_weights(p, fgeomap, fips, wips)

  # Face moments, i.e., M(Fi)_{ab} = w_Fi^b q_RF^a(xgp_RFi^b) (n_Fi × ())
  fshfs = QGradMonomialBasis{num_dims(fp)}(et,order-1)

  fmoments = Gridap.ReferenceFEs._Nedelec_face_moments(p, fshfs, c_fips, fcips, fwips)

  return fcips, fmoments

end

function Gridap.ReferenceFEs._Nedelec_cell_values(p,et,order)
  println("Cells: increased quad")

  # Compute integration points at interior
  degree = 2*(order) + 2
  iquad = Quadrature(p,degree)
  ccips = get_coordinates(iquad)
  cwips = get_weights(iquad)

  # Cell moments, i.e., M(C)_{ab} = q_C^a(xgp_C^b) w_C^b ⋅ ()
  if is_n_cube(p)
    cbasis = QCurlGradMonomialBasis{num_dims(p)}(et,order-1)
  else
    D = num_dims(p)
    cbasis = MonomialBasis{D}(VectorValue{D,et},order-D+1,(e,k)->sum(e)<=k)
  end
  cmoments = Gridap.ReferenceFEs._Nedelec_cell_moments(p, cbasis, ccips, cwips )

  return [ccips], [cmoments]

end



function uX(p)
  function _u(γαβ)
    xyz = forward_map_3D(p)(γαβ)
    # VectorValue(xyz[1]*xyz[2]*xyz[3], 0.0, 0.0)

    r = sqrt(xyz[1]^2 + xyz[2]^2 + xyz[3]^2)
    u = -xyz[2]/r
    v = xyz[1]/r
    w = 0.0
    VectorValue(u,v,w)
  end
end


# panel_model = models[1]
# p_fe = 0
function hodge_laplacian(
  panel_model::GridapDistributed.GenericDistributedDiscreteModel{3,3},
  p_fe::Int,dir::String,uX::Function,ls=LUSolver(),return_vtk=false)

  ranks = get_ranks(panel_model)

  lvl_h = nref(nc_horizontal(panel_model))
  lvl_v = nref(nc_vertical(panel_model))
  i_am_main(ranks) && println("nref_h = $lvl_h; nref_v = $lvl_v; p_fe = $p_fe")

  degree = 4*(p_fe + 1)
  if p_fe == 0
    degree = 8
  end
  @check degree > 0 "Zero quad!!"

  ## finite element solver
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*degree)
  Ω_error = Triangulation(panel_model)
  dΩ_error = Measure(Ω_error,4*degree)


  ## metric information
  inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  # manufacture rhs
  W(p) = x -> transpose_jacobian(p)(x) ⋅ surfcurl(uX,p)(x) # covariant comps of curl
  ccurl(p) = x-> 1/sqrtg(p,x)  * (forward_jacobian_3D(p,x) ⋅ curl(W(p))(x))
  cov_ccurl(p) = x -> transpose_jacobian(p)(x) ⋅ ccurl(p)(x)

  ccurl_cf = panelwise_cellfield(ccurl,Ω_panel,panel_ids)
  cov_ccurl_cf = panelwise_cellfield(cov_ccurl,Ω_panel,panel_ids)

  q_sdiv(p) = αβ ->  sqrtg(p,αβ)*( contra_v_3D(uX,p)(αβ))
  q_surfdiv(p::Int) = αβ -> 1/sqrtg(p,αβ) * ( divergence(q_sdiv(p))(αβ) )
  graddiv(p) = x -> sgrad(q_surfdiv,p)(x)
  cov_graddiv(p) = x -> transpose_jacobian(p)(x) ⋅ graddiv(p)(x)

  sdiv_cf = panelwise_cellfield(q_surfdiv,Ω_panel,panel_ids)

  rhs(p) = x -> ccurl(p)(x) - graddiv(p)(x)
  cov_rhs(p) = x -> cov_ccurl(p)(x) -  cov_graddiv(p)(x)

  rhs_cf = panelwise_cellfield(rhs,Ω_panel,panel_ids)
  cov_rhs_cf = panelwise_cellfield(cov_rhs,Ω_panel,panel_ids)

  graddiv_cf = panelwise_cellfield(cov_graddiv,Ω_panel,panel_ids)

  writevtk(Ω_panel,dir*"/rhs",
            cellfields=["w1"=>rhs_cf,"w2"=>cov_rhs_cf],
            append=false,geo_map= geo_map_func(Ω_panel))

  ## cellfields
  u_cov_cf = panelwise_cellfield(covar_v_3D(uX),Ω_panel,panel_ids)
  # sdiv_cf = panelwise_cellfield(surfdiv(contra_v_3D(uX)),Ω_panel,panel_ids)
  sigma_cf = -sdiv_cf


  tags = ["top_boundary", "bottom_boundary"]

  ## FE spaces
  T = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1, constraint=:zeromean)
  S = TrialFESpace(T)

  R = TestFESpace(Ω_panel, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl,dirichlet_tags=tags)
  H = TrialFESpace(R,u_cov_cf)

  X = MultiFieldFESpace([S,H])
  Y = MultiFieldFESpace([T,R])


  biform_s((s,u),(t,v)) = ( ∫( (s*t)*meas_cf )dΩ
                          - ∫( ( ∇(t)⋅(inv_metric_cf⋅u) )*meas_cf)dΩ
                         )
  biform_u((s,u),(t,v)) = ( ∫( (∇×u)⋅(metric_cf ⋅(∇×v) )*(1/meas_cf)  )dΩ
                          + ∫( ∇(s)⋅(inv_metric_cf⋅v )*meas_cf )dΩ
                          )
  biformX((s,u),(t,v)) = biform_s((s,u),(t,v)) + biform_u((s,u),(t,v))
  liformX((t,v)) = ∫( cov_rhs_cf⋅(inv_metric_cf⋅v)*meas_cf  )dΩ
  op = AffineFEOperator(biformX,liformX,X,Y)

  # using LinearAlgebra
  # A = get_matrix(op)
  # evals = eigvals(Array(partition(A).item))
  # A1 = partition(A).item
  # norm(A1-A1')/norm(A1)




  # dX=get_trial_fe_basis(X)
  # dY=get_fe_basis(Y)

  # A11((s,u),(t,v)) = ∫( (s*t)*meas_cf )dΩ
  # _A11=assemble_matrix(A11(dX,dY), X, Y)
  # A1 = partition(_A11).item
  # norm(A1-A1')/norm(A1)

  # A22((s,u),(t,v)) = ∫( (∇×u)⋅(metric_cf ⋅(∇×v) )*(1/meas_cf)  )dΩ
  # _A22=assemble_matrix(A22(dX,dY), X, Y)
  # A2 = partition(_A22).item
  # norm(A2-A2')/norm(A2)


  # A12((s,u),(t,v)) =  ∫( -1.0*( ∇(t)⋅(inv_metric_cf⋅u) )*meas_cf)dΩ
  # _A12=assemble_matrix(A12(dX,dY), X, Y)
  # a12 = partition(_A12).item
  # norm(a12-a12')/norm(a12)

  # A21((s,u),(t,v)) = ∫( ∇(s)⋅(inv_metric_cf⋅v )*meas_cf )dΩ
  # _A21=assemble_matrix(A21(dX,dY), X, Y)
  # a21 = partition(_A21).item
  # norm(a21-a21')/norm(a21)

  # a = a12 + a21'
  # norm(a)



  # using LinearAlgebra
  # A = get_matrix(op)
  # evals = eigvals(Array(partition(A).item))

  # A1 = partition(A).item
  # norm(A1-A1')/norm(A1)

  sh, uh = solve(ls,op)




  u_cov_int = interpolate(u_cov_cf,H)
  uh_ambient = covarient_basis_cf ⋅ (inv_metric_cf ⋅ uh )
  u_ambient = covarient_basis_cf ⋅ (inv_metric_cf ⋅ u_cov_int )

  _e = ( inv_metric_cf ⋅ uh) - ( inv_metric_cf ⋅  u_cov_int)
  el2_u =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ_error))

  _e = sh - sigma_cf
  el2_s =  sqrt(sum(∫(  (_e*_e)*meas_cf )dΩ_error))


  if return_vtk

    cellfields = ["uambient"=>u_ambient,
                  "uambient_h"=>uh_ambient,
                  "eu_ambient"=>u_ambient-uh_ambient,
                  "eu"=>u_cov_cf-uh,
                  "s"=>sigma_cf,
                  "sh"=>sh,
                  "es"=>sigma_cf-sh  ]

    ### plot in 3D
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl_h)_p$p_fe",
            cellfields=cellfields,append=false,geo_map= geo_map_func(Ω_panel))


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



  return el2_u, false, false
  # return el2_u, el2_s, false

end


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

dir = datadir("HodgeLaplacianConvergence")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

n_ref_lvls = 4
ps = [0,1]
models  = get_3D_octree_refined_models(ranks,n_ref_lvls)
ls = LUSolver()
p_convergence_test(ranks,ps,models,hodge_laplacian,dir,uX,ls,true)
# i_am_main(ranks) && plot_convergence_from_saved(dir,"convergence_p0",["u:", "s:"])
