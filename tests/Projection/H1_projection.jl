using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed
using DrWatson

include("../convergence_tools.jl")


### f = sin(ϕ) = Z/R
function func(p)
  function _f(αβ)
    X,Y,Z = ForwardMap(p)(αβ)
    radius = sqrt(X^2 + Y^2 + Z^2)
    Z/radius
  end
end


MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))


n_ref_lvls = 4
ps = [1,2]

dir = datadir("InterpolationConvergence_glued")
!isdir(dir) && mkdir(dir)


models = get_octree_refined_models(ranks,n_ref_lvls)



function interpolation(panel_model,p_fe::Int,dir::String,func::Function,return_vtk)


  Dc = num_cell_dims(panel_model)

  lvl = nref(nc(panel_model))


  println("p_fe = $(p_fe); nref = $lvl; Dc = $Dc")

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,5*p_fe)
  panel_ids = get_panel_ids(panel_model)

  f_panel_cf = panelwise_cellfield(func,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)

  V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  U = TrialFESpace(V)


  a(u,v) = ∫( (u*v)*meas_cf )dΩ + ∫( ( ∇(v)⋅(inv_metric_cf⋅∇(u)))*meas_cf )dΩ
  l(v) = ∫( (f_panel_cf*v)*meas_cf )dΩ + ∫( ( ∇(v)⋅(inv_metric_cf⋅∇(f_panel_cf)))*meas_cf )dΩ
  op = AffineFEOperator(a,l,U,V)
  f_uh = solve(LUSolver(),op)


  _e = f_panel_cf - f_uh
  el2  =  sqrt( sum(∫( (_e*_e)*meas_cf )dΩ) )
  eh1  =  sqrt( sum( ∫( (_e*_e)*meas_cf )dΩ + ∫( (∇(_e)⋅ (inv_metric_cf⋅∇(_e)) )*meas_cf )dΩ   )  )

  if return_vtk
    panel_cfs = [f_panel_cf, f_uh,  _e, gradient(f_uh) ]
    labels = ["u","uh", "e" , "grad"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$(p_fe)_"*String(:H1),
            cellfields=cellfields,append=false,geo_map=latlon_geo_map_func(Ω_panel))
  end

  return eh1,el2,false
end


models = get_octree_refined_models(ranks,n_ref_lvls)
_dir = dir*"/H1_interpolate_scalar_func_2D"
!isdir(_dir) && mkdir(_dir)
p_convergence_test(ranks,ps,models,interpolation,_dir,func,true)
plot_convergence_from_saved(_dir,"convergence",["h1", "l2"])
