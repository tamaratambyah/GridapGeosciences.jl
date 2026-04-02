using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed
using DrWatson

include("../convergence_tools.jl")


function interpolation(panel_model,p_fe::Int,dir::String,func::Function,conf,return_vtk)

  Dc = num_cell_dims(panel_model)

  lvl = 0
  if Dc == 2
    lvl = nref(nc(panel_model))
  elseif Dc == 3
    lvl = nref(nc_horizontal(panel_model))
  end

  println("p_fe = $(p_fe); nref = $lvl; Dc = $Dc")

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,5*p_fe)
  panel_ids = get_panel_ids(panel_model)

  f_panel_cf = panelwise_cellfield(func,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)

  V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  U = TrialFESpace(V)

  if Dc == 3
    V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=conf,
                    dirichlet_tags=["top_boundary", "bottom_boundary"])
    U = TrialFESpace(V,f_panel_cf)
  end

  # f_uh = interpolate(f_panel_cf,U)
  a(u,v) = ∫( (u*v)*meas_cf )dΩ
  l(v) = ∫( (f_panel_cf*v)*meas_cf )dΩ
  op = AffineFEOperator(a,l,U,V)
  f_uh = solve(LUSolver(),op)

  _e = f_panel_cf - f_uh
  e  =  sqrt( sum(∫( (_e*_e)*meas_cf )dΩ) )

  if return_vtk
    panel_cfs = [f_panel_cf, f_uh,  _e, gradient(f_uh) ]
    labels = ["u","uh", "e" , "grad"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$(p_fe)_"*String(conf),
            cellfields=cellfields,append=false,geo_map=latlon_geo_map_func(Ω_panel))
  end

  return e,false,false

end

# function fS(p)
#   function f(αβ)
#     xyz = forward_map_2D(p)(αβ)
#     xyz[1]*xyz[2]*xyz[3]
#   end
# end

### f = sin(ϕ) = Z/R
function fS(p)
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

Dc = 2

models = (Dc == 2) ? get_octree_refined_models(ranks,n_ref_lvls) : get_3D_octree_refined_models(ranks,n_ref_lvls)
for conf in [:L2, :H1]
  _dir = dir*"/scalar_func_$(Dc)D_"*String(conf)
  !isdir(_dir) && mkdir(_dir)
  p_convergence_test(ranks,ps,models,interpolation,_dir,fS,conf,true)
  plot_convergence_from_saved(_dir,"convergence",["p"])
end
