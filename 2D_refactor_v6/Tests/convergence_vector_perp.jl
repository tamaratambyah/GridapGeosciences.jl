################################################################################
#### Test normal vectors
################################################################################
function test_normal_unit_vector(panel_model,return_vtk=false)
  Ω_panel = Triangulation(panel_model)
  panel_ids = get_panel_ids(panel_model)
  dΩ = Measure(Ω_panel,4)

  vX = panel_to_cartesian(normal_vec)
  norm_vec_cf = panelwise_cellfield(vX,Ω_panel,panel_ids)
  norm_vec_from_basis_cf = panelwise_cellfield(normal_vector_from_basis,Ω_panel,panel_ids)

  @check l2(norm_vec_cf-norm_vec_from_basis_cf,dΩ) < 1e-16

  if return_vtk
    lvl = nref(nc(panel_model))
    cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)
    panel_cfs = [ norm_vec_cf,norm_vec_from_basis_cf]
    labels = ["normal", "n"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$lvl",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end
end

models  = get_refined_models(n_ref_lvls)
for panel_model in models
  test_normal_unit_vector(panel_model)
end


################################################################################
#### Perp convergence
################################################################################

function vector_perp(panel_model,vecX::Function,p_fe::Int,return_vtk=false)
  Ω_panel = Triangulation(panel_model)
  panel_ids = get_panel_ids(panel_model)
  dΩ = Measure(Ω_panel,4*p_fe)


  vX = panel_to_cartesian(tangent_vec(vecX))

  norm_vec_from_basis_cf = panelwise_cellfield(normal_vector_from_basis,Ω_panel,panel_ids)
  u_proj_cf = panelwise_cellfield(projection_v(vX),Ω_panel,panel_ids)
  u_perp = cross(norm_vec_from_basis_cf,u_proj_cf)

  meas_cf = CellField(sqrtg,Ω_panel)
  jacobian_cf = panelwise_cellfield(forward_jacobian,Ω_panel,panel_ids)
  u_perp_contra = panelwise_cellfield(contra_v_perp(vX),Ω_panel,panel_ids)

  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  u_perp_contrah = interpolate(u_perp_contra,U)
  u_perph = jacobian_cf ⋅ u_perp_contrah

  e = l2((u_perp - u_perph)*meas_cf,dΩ)

  if return_vtk
    lvl = nref(nc(panel_model))
    cell_geo_map = lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), panel_ids)

    panel_cfs = [u_perph, u_perp,u_proj_cf,u_perph-u_perp ]
    labels = ["u_perph","u_perp","u_proj", "e"]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$lvl",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  return e

end



function vector_perp_errors(panel_model,func::Function,p_fe::Int,return_vtk=false)
  e  = vector_perp(panel_model,func,p_fe,return_vtk)
  return e,false
end

function vector_perp_convergence_test(analytic_funcs,n_ref_lvls,return_vtk=false)

  for (key, val) in analytic_funcs
    plot()
    for p_fe in [1,2,3]
      errs,ns,dxs,slope = convergence_test(vector_perp_errors,n_ref_lvls,val,p_fe,return_vtk)
      plot_convergence(errs,ns,dxs,slope;leginf=["p=$p_fe"],colors=[palette(:tab10)[p_fe]])
    end
    savefig(plotsdir()*"/vector_perp_convergence_func_$(key)")
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

vector_perp_convergence_test(ambient_vecs,n_ref_lvls,false)



### latlon vectors
vecθϕ_1(θϕ) = VectorValue(cos(θϕ[1]),0.0)
vecθϕ_2(θϕ) = VectorValue(-sin(θϕ[1]),0.0)

latlon_vecs = Dict{Symbol,Any}()
latlon_vecs[:vtheta1] = vec_cartesian_to_latlon(vecθϕ_1)
# latlon_vecs[:vtheta2] = vec_cartesian_to_latlon(vecθϕ_2)
n_ref_lvls = 4

vector_perp_convergence_test(latlon_vecs,n_ref_lvls,true)


### williamson2 vector field
vWilliamson(ζ) = θϕ -> - VectorValue( cos(θϕ[2])*cos(ζ) + cos(θϕ[1])*sin(θϕ[2])*sin(ζ),
                                      -sin(θϕ[1])*sin(ζ) )
williamson_vec = Dict{Symbol,Any}()
# williamson_vec[:z1] = vec_cartesian_to_latlon(vWilliamson(0.0))
# williamson_vec[:z2] = vec_cartesian_to_latlon(vWilliamson(0.05))
# williamson_vec[:z3] = vec_cartesian_to_latlon(vWilliamson(π/2-0.05))
williamson_vec[:z4] = vec_cartesian_to_latlon(vWilliamson(π/2))
vector_perp_convergence_test(williamson_vec,n_ref_lvls,true)
