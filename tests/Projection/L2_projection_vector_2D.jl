using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed
using DrWatson

include("../convergence_tools.jl")


function interpolation(panel_model::GridapDistributed.GenericDistributedDiscreteModel{2,2},
                        p_fe::Int,dir::String,vecX::Function,conf,return_vtk)
  Dc = num_cell_dims(panel_model)

  lvl = nref(nc(panel_model))
  println("p_fe = $(p_fe); nref = $lvl; Dc = $Dc")

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,4*p_fe)
  dΩ_error = Measure(Ω_panel,8*p_fe)
  panel_ids = get_panel_ids(panel_model)

  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  vec_contra_cf = panelwise_cellfield(contra_v(vecX),Ω_panel,panel_ids)
  vec_proj_cf = covarient_basis_cf⋅vec_contra_cf

  reffe = ReferenceFE(lagrangian,VectorValue{Dc,Float64},p_fe)
  if conf == :HDiv || conf == :Hdiv
    reffe = ReferenceFE(raviart_thomas,Float64,p_fe)
  end

  V = TestFESpace(panel_model, reffe; conformity=conf)
  U = TrialFESpace(V)
  # vec_contra_h = interpolate(vec_contra_cf,U)

  a(u,v) = ∫( u⋅( metric_cf⋅v)*meas_cf )dΩ
  l(v) = ∫( vec_contra_cf⋅(metric_cf⋅v)*meas_cf )dΩ
  op = AffineFEOperator(a,l,U,V)
  vec_contra_h = solve(LUSolver(),op)

  vec_proj_h = covarient_basis_cf ⋅vec_contra_h

  _e = vec_contra_cf - vec_contra_h
  el2_proj =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ_error))


  vec_contra_h  = interpolate(vec_contra_cf, V)
  vec_interp_h = covarient_basis_cf ⋅vec_contra_h

  _e = vec_contra_cf - vec_contra_h
  el2_interp =  sqrt(sum(∫( _e⋅(metric_cf⋅_e)*meas_cf )dΩ_error))



  if return_vtk
    panel_cfs = [vec_proj_cf, vec_proj_h, vec_proj_cf-vec_proj_h,
                vec_interp_h, vec_interp_h-vec_proj_cf ]
    labels = ["u_proj", "u_projh", "eproj",
            " u_int", "e_int"]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$(p_fe)_"*String(conf),cellfields=cellfields,
          append=false,geo_map=latlon_geo_map_func(Ω_panel))
  end

  return el2_interp,el2_proj,false
end

### must be in the tangent space of the sphere
vX(xyz) = VectorValue(-xyz[2], xyz[1], 0)
vecX = panel_to_cartesian(tangent_vec(vX))



MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))


n_ref_lvls = 4
ps = [2]

dir = datadir("InterpolationConvergence")
!isdir(dir) && mkdir(dir)

Dc = 2
models  = get_octree_refined_models(ranks,n_ref_lvls)

for conf in [:Hdiv]
  _dir = dir*"/vector_func_$(Dc)D_"*String(conf)
  !isdir(_dir) && mkdir(_dir)
  p_convergence_test(ranks,ps,models,interpolation,_dir,vecX,conf,true)
  plot_convergence_from_saved(_dir,"convergence",["Interp","L2Proj", ])
end
