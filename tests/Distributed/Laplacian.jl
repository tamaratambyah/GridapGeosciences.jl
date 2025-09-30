
################################################################################
#### DistributedParametricDiscreteModel
################################################################################

a_dmodels = Vector{DistributedParametricDiscreteModel}(undef,length(models))
dpanel_ids = Vector{AbstractArray{Vector{Int}}}(undef,length(models))
owned_panel_ids = Vector{AbstractArray{Vector{Int}}}(undef,length(models))

coarse_dmodel = DiscreteModel(ranks,models[end],cell_to_part[end])
dpanel_ids[end],owned_panel_ids[end] = distributed_panel_ids(coarse_dmodel,spanel_ids[end])

a_dmodels[end] = DistributedParametricDiscreteModel(coarse_dmodel,dpanel_ids[end])
for level in length(models)-1:-1:1
  child = DiscreteModel(ranks,models[level],cell_to_part[level])
  parent = a_dmodels[level+1]
  glue = DistributedAdaptivityGlue(glues[level],parent,child)
  dpanel_ids[level],owned_panel_ids[level] = distributed_panel_ids(child,spanel_ids[level])

  a_dmodels[level] = DistributedAdaptedParametricDiscreteModel(child,parent,glue,dpanel_ids[level])
end


get_panel_ids(a_dmodels[1])

get_model(a_dmodels[1])


################################################################################
#### Laplacian
################################################################################
include("../Laplace/LaplaceBeltrami.jl")
include("../Laplace/analytic_funcs.jl")

dir = datadir("DistributedLaplaceTests")
return_vtk = true
(return_vtk && !isdir(dir)) && mkdir(dir)

n_ref_lvls = 4

nprocs = 6
ranks  = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end

coarse_s_model = true
dmodels, dpanel_ids, owned_panel_ids = get_distributed_refined_models(ranks,nprocs,n_ref_lvls,coarse_s_model)

p_fe = 1
f =  f_XYZ
panel_model = dmodels[1]
e, uh, f_panel_cf = laplace_beltrami_solver(panel_model,f,p_fe,return_vtk)



function laplace_beltrami_convergence_test(ranks,analytic_funcs,n_ref_lvls,return_vtk=false)

  for (key, val) in analytic_funcs
    for p_fe in [1, 2, 3]
      errs,ns,dxs,slope = convergence_test(laplace_beltrami_errors,n_ref_lvls,val,p_fe,return_vtk)
      output = @strdict errs ns dxs slope
      safesave(datadir(dir, ("laplace_beltrami_convergence_func_$(key)_p$(p_fe).jld2")), output)
    end
  end

end



cell_geo_map = geo_map_func(panel_ids)

# cell_geo_map = map(panel_ids) do pid
#   return lazy_map(p -> MatMultField(R1p[p]) ∘ ForwardMapPanel1(), pid)
# end


writevtk(Triangulation(panel_model),dir*"/ambient_model_cf_ref_$(level)",cellfields=["f"=>f_panel_cf],
append=false,geo_map=cell_geo_map)



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
eh = f_panel_cf-uh

writevtk(Ω_panel,dir*"/ambient_model_nref$(level)_p$p_fe",cellfields=["eu"=>eh],append=false,geo_map=cell_geo_map)
