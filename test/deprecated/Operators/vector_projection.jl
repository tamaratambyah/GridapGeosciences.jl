"""
test the projection of ambient vector fields on the tangent space
"""


function vector_proj(panel_model,p_fe::Int,dir::String,vecX::Function,return_vtk=false)
  lvl = nref(nc(panel_model))
  println("p_fe = $(p_fe); nref = $lvl")

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,4*p_fe)
  panel_ids = get_panel_ids(panel_model)

  vec_phys = ParametricCellField(vecX,Ω_panel,panel_ids)
  vec_project = ParametricCellField(projection_v(vecX),Ω_panel,panel_ids)
  vec_contra_cf = ParametricCellField(contra_v(vecX),Ω_panel,panel_ids)
  jacobian_cf = ParametricCellField(forward_jacobian,Ω_panel,panel_ids)

  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)
  vec_contra_h = interpolate(vec_contra_cf,U)
  project_h = jacobian_cf ⋅vec_contra_h

  e = l2(project_h - vec_project,dΩ)

  if return_vtk
    lvl = nref(nc(panel_model))
    cell_geo_map = geo_map_func(Ω_panel)

    panel_cfs = [vec_phys, vec_project, project_h, project_h - vec_project ]
    labels = ["u","u_proj", "u_projh", "e"]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  return e, false,false

end


function vector_proj_convergence_test(ranks::AbstractArray,nprocs::Int,analytic_funcs,
    n_ref_lvls=4,ps=[1,2,3],return_vtk=false)
  # serial models
  models  = get_refined_models(n_ref_lvls)
  dir = datadir("VectorProjectionConvergence")
  !isdir(dir) && mkdir(dir)

  for (key, val) in analytic_funcs
    _dir = dir*"/func_$(key)"
    !isdir(_dir) && mkdir(_dir)
    p_convergence_test(ranks,ps,models,vector_proj,_dir,val,return_vtk)
    plot_convergence_from_saved(_dir,"convergence",["p"])
  end

end
