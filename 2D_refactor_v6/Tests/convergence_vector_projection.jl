
function vector_proj(panel_model,vecX::Function,p_fe::Int,return_vtk=false)

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,4*p_fe)
  panel_ids = get_panel_ids(panel_model)

  vec_phys = panelwise_cellfield(vecX,Ω_panel,panel_ids)
  vec_project = panelwise_cellfield(projection_v(vecX),Ω_panel,panel_ids)
  vec_contra_cf = panelwise_cellfield(contra_v(vecX),Ω_panel,panel_ids)
  jacobian_cf = panelwise_cellfield(forward_jacobian,Ω_panel,panel_ids)

  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)
  vec_contra_h = interpolate(vec_contra_cf,U)
  project_h = jacobian_cf ⋅vec_contra_h

  e = l2(project_h - vec_project,dΩ)

  if return_vtk
    lvl = nref(nc(panel_model))
    cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)

    panel_cfs = [vec_phys, vec_project, project_h, project_h - vec_project ]
    labels = ["u","u_proj", "u_projh", "e"]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$lvl",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  return e, project_h, vec_project

end



function vector_proj_errors(panel_model,func::Function,p_fe::Int,return_vtk=false)
  e,  = vector_proj(panel_model,func,p_fe,return_vtk)
  return e,false
end

function vector_proj_convergence_test(analytic_funcs,n_ref_lvls,return_vtk=false)

  for (key, val) in analytic_funcs
    plot()
    vX = panel_to_cartesian(tangent_vec(val)) # ensure tangent vector
    for p_fe in [1,2,3]
      errs,ns,dxs,slope = convergence_test(vector_proj_errors,n_ref_lvls,vX,p_fe,return_vtk)
      plot_convergence(errs,ns,dxs,slope;leginf=["p=$p_fe"],colors=[palette(:tab10)[p_fe]])
    end
    savefig(plotsdir()*"/vector_proj_convergence_func_$(key)")
  end

end


### ambient vectors
vecX_1(XYZ) = VectorValue(-1.0*XYZ[2],XYZ[1],0.0)
vecX_2(XYZ) = VectorValue(XYZ[1]*XYZ[2],XYZ[2]*XYZ[3],XYZ[3]^2-RADIUS^2)
vecX_3(XYZ) = VectorValue(XYZ[2],XYZ[3],0.0)

ambient_vecs = Dict{Symbol,Any}()
ambient_vecs[:v1] = vecX_1
ambient_vecs[:v2] = vecX_2
ambient_vecs[:v3] = vecX_3
n_ref_lvls = 4

vector_proj_convergence_test(ambient_vecs,n_ref_lvls,false)


### latlon vectors
vecθϕ_1(θϕ) = VectorValue(cos(θϕ[1]),0.0)
vecθϕ_2(θϕ) = VectorValue(-sin(θϕ[1]),0.0)

latlon_vecs = Dict{Symbol,Any}()
# latlon_vecs[:vtheta1] = vec_cartesian_to_latlon(vecθϕ_1)
latlon_vecs[:vtheta2] = vec_cartesian_to_latlon(vecθϕ_2)
n_ref_lvls = 4

vector_proj_convergence_test(latlon_vecs,n_ref_lvls,true)


### williamson2 vector field
vWilliamson(ζ) = θϕ -> - VectorValue( cos(θϕ[2])*cos(ζ) + cos(θϕ[1])*sin(θϕ[2])*sin(ζ),
                                      -sin(θϕ[1])*sin(ζ) )
williamson_vec = Dict{Symbol,Any}()
williamson_vec[:z1] = vec_cartesian_to_latlon(vWilliamson(0.0))
# williamson_vec[:z2] = vec_cartesian_to_latlon(vWilliamson(0.05))
# williamson_vec[:z3] = vec_cartesian_to_latlon(vWilliamson(π/2-0.05))
# williamson_vec[:z4] = vec_cartesian_to_latlon(vWilliamson(π/2))
vector_proj_convergence_test(williamson_vec,n_ref_lvls,true)
