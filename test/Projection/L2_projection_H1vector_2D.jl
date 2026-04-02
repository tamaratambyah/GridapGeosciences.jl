using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed
using DrWatson
using LinearAlgebra
using FillArrays

include("../convergence_tools.jl")


function interpolation(
  dpanel_model::Union{GridapDistributed.GenericDistributedDiscreteModel{2,2},GridapDistributed.GenericDistributedDiscreteModel{3,3}},
                        p_fe::Int,dir::String,vecX::Function,return_vtk)

  panel_model = dpanel_model
  ranks = get_ranks(panel_model)

  Dc = num_cell_dims(panel_model)

  lvl = 0
  if Dc == 2
    lvl = nref(nc(panel_model))
  elseif Dc == 3
    lvl = nref(nc_horizontal(panel_model))
  end

  i_am_main(ranks) && println("p_fe = $(p_fe); nref = $lvl; Dc = $Dc")

  degree = 4*(p_fe+1)

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)
  dΩ_error = Measure(Ω_panel,2*degree)
  panel_ids = get_panel_ids(panel_model)

  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  vec_contra_cf = panelwise_cellfield(contra_v(vecX),Ω_panel,panel_ids)
  vec_proj_cf = covarient_basis_cf⋅vec_contra_cf


  a(u,v) = ∫( (u⋅(metric_cf⋅v))*meas_cf )dΩ
  l(v) = ∫( (vec_contra_cf⋅(metric_cf⋅v))*meas_cf )dΩ

  reffe  = ReferenceFE(lagrangian,VectorValue{Dc, Float64},p_fe)
  V = TestFESpace(panel_model, reffe; conformity=:H1)
  U = TrialFESpace(V)

  if Dc == 3
    V = TestFESpace(panel_model, reffe; conformity=:H1,
                dirichlet_tags=["top_boundary", "bottom_boundary"])
    U = TrialFESpace(V,vec_contra_cf)
  end

  op = AffineFEOperator(a,l,U,V)

  vec_contra_h = solve(op)
  vec_l2proj_h = covarient_basis_cf ⋅vec_contra_h

  _e = vec_contra_cf - vec_contra_h
  el2_proj =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ_error))

  # Interpolation instead of L2 projection ...
  vec_contra_h = interpolate(vec_contra_cf, U)
  vec_interp_h = covarient_basis_cf ⋅vec_contra_h
  _e = vec_contra_cf - vec_contra_h
  el2_interp =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ_error))

  i_am_main(ranks) && println("Error interp: ", el2_interp)
  i_am_main(ranks) && println("Error proj: ", el2_proj)

  if return_vtk
    panel_cfs = [vec_proj_cf, vec_l2proj_h, vec_proj_cf-vec_l2proj_h,
                vec_interp_h, vec_interp_h-vec_proj_cf]
    labels = ["u_proj", "u_projh", "eproj",
              "u_int", "e_int"]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$(p_fe)",cellfields=cellfields,
          append=false,geo_map=geo_map_func(Ω_panel))
  end

  return  el2_interp,el2_proj,false

end

### must be in the tangent space of the sphere
function uX(p)
  function _u(α)
    x = ForwardMap(p)(α)
    VectorValue(-x[2],x[1],0.0)
  end
end

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

n_ref_lvls = 4
ps = [1]

dir = datadir("InterpolationConvergence")
!isdir(dir) && mkdir(dir)

Dc = 3
models = (Dc == 2) ? get_octree_refined_models(ranks,n_ref_lvls) : get_3D_octree_refined_models(ranks,n_ref_lvls-1)

_dir = dir*"/vector_func_$(Dc)D_H1"
!isdir(_dir) && mkdir(_dir)
p_convergence_test(ranks,ps,models,interpolation,_dir,uX,true)
plot_convergence_from_saved(_dir,"convergence",["Interp","L2Proj", ])
