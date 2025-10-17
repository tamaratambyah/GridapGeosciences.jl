""" Helmholtz problem in mixed form
σ - ∇ᵧ(u) = 0
u + ∇ᵧ⋅σ = f
"""

module MixedHelmholtz

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using GridapGeosciences
using Test

include("analytic_funcs.jl")
include("../convergence_tools.jl")

function mixed_helmholtz_solver(panel_model,p_fe::Int,dir::String,f::Function,ls=LUSolver(),return_vtk=false)
  ranks = get_ranks(panel_model)

  lvl = nref(nc(panel_model))
  i_am_main(ranks) && println("nref = $lvl")

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  f_panel_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  sigma_cf = panelwise_cellfield(contr_gradf(f),Ω_panel,panel_ids)
  sdiv_cf =  panelwise_cellfield(surfdiv(contr_gradf(f)),Ω_panel,panel_ids)
  slap_panel_cf =  panelwise_cellfield(surflap(f),Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  gradu_cf = covarient_basis_cf ⋅ sigma_cf

  rhs_cf = f_panel_cf + sdiv_cf

  i_am_main(ranks) && println("Check sdiv(sgrad) against slap: ", l2(sdiv_cf-slap_panel_cf,dΩ) )


  V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  U = TrialFESpace(V)

  T = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  S = TrialFESpace(T)

  Y = MultiFieldFESpace([V, T])
  X = MultiFieldFESpace([U, S])

  metric_cf = CellField(analytic_metric,Ω_panel)
  meas_cf = CellField(sqrtg,Ω_panel)
  grad_meas_cf = CellField(grad_meas,Ω_panel)

  biform1((u,s),(v,t)) = ∫( (s⋅ (metric_cf⋅t))*meas_cf )dΩ + ∫( u*(t⋅grad_meas_cf + meas_cf*(∇⋅t) ) )dΩ
  biform2((u,s),(v,t)) = ∫( (u*v)*meas_cf )dΩ + ∫( v*(s⋅grad_meas_cf + meas_cf*(∇⋅s) ) )dΩ

  biformX((u,s),(v,t)) = biform1((u,s),(v,t)) + biform2((u,s),(v,t))
  liformX((v,t)) = ∫( (rhs_cf*v)*meas_cf )dΩ

  op = AffineFEOperator(biformX,liformX,X,Y)
  # uh,sh = solve(ls,op)

  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = Gridap.Algebra.allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  xh = FEFunction(X,x)
  uh,sh = xh

  graduh = covarient_basis_cf ⋅sh

  e_u = l2( (f_panel_cf - uh)*meas_cf,dΩ) # error in scalar u
  e_s = l2((sigma_cf - sh)*meas_cf,dΩ) # error in contra compons of grad u
  e_gradu = l2((gradu_cf - graduh)*meas_cf,dΩ) # error in grad u = physical sigma

  if return_vtk
    lvl = nref(nc(panel_model))
    cell_geo_map = geo_map_func(Ω_panel)

    panel_cfs = [uh, sh, graduh, f_panel_cf, gradu_cf, rhs_cf, f_panel_cf - uh, sigma_cf - sh, gradu_cf - graduh  ]
    labels = ["uh","sh","graduh", "u_ex", "gradu","rhs", "eu", "es", "egradu"]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  return e_u,e_gradu,false

end

################################################################################
#### Auto convergence test
################################################################################
function main(distribute,nprocs;octree=false)
  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("Auto conference test: Mixed Helmholtz")

  n_ref_lvls = 4
  ps = [1,2] ## using RT, so only test p=1,2
  ls = LUSolver()
  models  = get_refined_models(n_ref_lvls)

  dir = datadir("MixedHelmholtzConvergence")
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  if prod(nprocs) > 1
    i_am_main(ranks) && println("Distributed test")
    if octree
      i_am_main(ranks) && println("Octrees")
      models =  get_octree_refined_models(ranks,n_ref_lvls)
    else
      models,  = get_distributed_refined_models(ranks,nprocs,models)
    end
    ## can use minres solver, but it takes ages to converge. So just use LU
    # ls = MINRESSolver(;Pl=JacobiLinearSolver(),maxiter=5000,verbose=true)
  end

  for (key, val) in analytic_funcs
    i_am_main(ranks) && println("mixed_helmholtz_convergence_func_$(key)")
    p_convergence_test(ranks,ps,models,mixed_helmholtz_solver,dir,val,ls)
  end

  i_am_main(ranks) && println("--DONE--")
end

################################################################################
#### Convergence test with plots
################################################################################

function mixed_helmholtz_convergence_test(ranks::AbstractArray,nprocs::Int,
  analytic_funcs,n_ref_lvls=4,ps=[1],ls=LUSolver(),return_vtk=false)

  # serial models
  models  = get_refined_models(n_ref_lvls)

  if prod(nprocs) > 1
    i_am_main(ranks) && println("Distributed test")
    models,  = get_distributed_refined_models(ranks,nprocs,models)
  end

  dir = datadir("MixedHelmholtz")
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  for (key, val) in analytic_funcs
    simName = "mixed_helmholtz_convergence_func_$(key)"
    i_am_main(ranks) && println(simName)

    errors = Vector{Vector{Float64}}(undef,length(ps))
    ns = Vector{Vector{Float64}}(undef,length(ps))
    dxs = Vector{Vector{Float64}}(undef,length(ps))
    slopes = Vector{Float64}(undef,length(ps))


    for (i,p_fe) in enumerate(ps)
      i_am_main(ranks) && println("p_fe = $p_fe")
      errors[i],ns[i],dxs[i],slopes[i] = h_convergence_test(models,mixed_helmholtz_solver,p_fe,dir,val,ls,return_vtk)
    end
    print_convergence_results(errors,ns,dxs,slopes,ps)
    output = @strdict errors ns dxs slopes ps

    i_am_main(ranks) && safesave(datadir(dir, ("$simName.jld2")), output)

    i_am_main(ranks) && plot_convergence_from_saved(dir,simName,["u","s"])

  end

end



end
