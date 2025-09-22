""" Poisson problem using Laplace-Beltrami operator
u + Δᵧ(u) = f
Need to remove the kernal via zeromean FE space
"""

using Gridap.Helpers
function laplace_beltrami_solver(panel_model,f::Function,p_fe::Int,return_vtk=false)
  lvl = nref(nc(panel_model))
  println("nref = $lvl")

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*p_fe+1)

  V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1, constraint=:zeromean)
  U = TrialFESpace(V)

  f_panel_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  inv_metric_cf = CellField(analytic_inv_metric,Ω_panel)
  meas_cf = CellField(sqrtg,Ω_panel)
  slap_panel_cf =  panelwise_cellfield(surflap(f),Ω_panel,panel_ids)

  println("Zeromean: ", sum(∫(f_panel_cf*meas_cf)dΩ))
  @check sum(∫(f_panel_cf*meas_cf)dΩ) < 1e-14 "Function must be zero mean to solve with zeromean FE space!"

  rhs_cf = - slap_panel_cf

  poisson_biform(u,v) =  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_cf )dΩ
  poisson_liform(v) = ∫(  (rhs_cf*v)*meas_cf )dΩ
  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)


  uh = solve(LUSolver(),op)

  e = l2(f_panel_cf-uh,dΩ)

  if return_vtk
    lvl = nref(nc(panel_model))
    cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
    panel_cfs = [f_panel_cf,uh,f_panel_cf-uh]
    labels = ["u","uh","eu"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end
  return e, uh, f_panel_cf
end


function laplace_beltrami_errors(panel_model,func::Function,p_fe::Int,return_vtk=false)
  e,  = laplace_beltrami_solver(panel_model,func,p_fe,return_vtk)
  return e,false,false
end

function laplace_beltrami_convergence_test(analytic_funcs,n_ref_lvls,return_vtk=false)

  for (key, val) in analytic_funcs
    plot()
    for p_fe in [1, 2, 3]
      errs,ns,dxs,slope = convergence_test(laplace_beltrami_errors,n_ref_lvls,val,p_fe,return_vtk)
      plot_convergence(errs,ns,dxs,slope;leginf=["p=$p_fe"],colors=[palette(:tab10)[p_fe]])
    end
    savefig(plotsdir()*"/laplace_beltrami_convergence_func_$(key)")
  end

end
