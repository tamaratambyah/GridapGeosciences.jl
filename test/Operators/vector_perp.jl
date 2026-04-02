################################################################################
#### Perp convergence
################################################################################


function vector_perp(panel_model,p_fe::Int,dir::String,vX::Function,return_vtk=false)
  lvl = nref(nc(panel_model))
  println("p_fe = $(p_fe); nref = $lvl")

  Ω_panel = Triangulation(panel_model)
  panel_ids = get_panel_ids(panel_model)
  dΩ = Measure(Ω_panel,4*p_fe)

  norm_vec_from_basis_cf = panelwise_cellfield(normal_vector_from_basis,Ω_panel,panel_ids)
  u_proj_cf = panelwise_cellfield(projection_v(vX),Ω_panel,panel_ids)
  u_perp = cross(norm_vec_from_basis_cf,u_proj_cf)

  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  jacobian_cf = panelwise_cellfield(forward_jacobian,Ω_panel,panel_ids)
  u_perp_contra = panelwise_cellfield(contra_v_perp(vX),Ω_panel,panel_ids)

  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  u_perp_contrah = interpolate(u_perp_contra,U)
  u_perph = jacobian_cf ⋅ u_perp_contrah

  e = l2((u_perp - u_perph),meas_cf,dΩ)

  if return_vtk
    lvl = nref(nc(panel_model))
    cell_geo_map = geo_map_func(Ω_panel)

    panel_cfs = [u_perph, u_perp,u_proj_cf,u_perph-u_perp ]
    labels = ["u_perph","u_perp","u_proj", "e"]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  return e,false,false

end


function vector_perp_convergence_test(ranks::AbstractArray,nprocs::Int,analytic_funcs,
  n_ref_lvls=4,ps=[1,2,3],return_vtk=false)
  # serial models
  models  = get_refined_models(n_ref_lvls)
  dir = datadir("VectorPerpConvergence")
  !isdir(dir) && mkdir(dir)

  for (key, val) in analytic_funcs
    _dir = dir*"/func_$(key)"
    !isdir(_dir) && mkdir(_dir)
    p_convergence_test(ranks,ps,models,vector_perp,_dir,val,return_vtk)
    plot_convergence_from_saved(_dir,"convergence",["p"])
  end

end
