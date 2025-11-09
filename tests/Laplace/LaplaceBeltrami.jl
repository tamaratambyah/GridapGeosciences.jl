""" Poisson problem using Laplace-Beltrami operator
u + Δᵧ(u) = f
Need to remove the kernal via zeromean FE space
"""

module LaplaceBeltrami

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

function laplace_beltrami_solver(
  panel_model::Union{<:DiscreteModel{2,2},<:GridapDistributed.DistributedDiscreteModel{2,2}},
  p_fe::Int,dir::String,f::Function,ls=LUSolver(),return_vtk=false)

  ranks = get_ranks(panel_model)

  lvl = nref(nc(panel_model))
  i_am_main(ranks) && println("nref = $lvl; p_fe = $p_fe")

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*p_fe+1)

  V = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1, constraint=:zeromean)
  U = TrialFESpace(V)

  f_panel_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  slap_panel_cf =  panelwise_cellfield(surflap(f),Ω_panel,panel_ids)

  # i_am_main(ranks) && println("Zeromean: ", sum(∫(f_panel_cf*meas_cf)dΩ))
  @check sum(∫(f_panel_cf*meas_cf)dΩ) < 1e-14 "Function must be zero mean to solve with zeromean FE space!"

  rhs_cf = - slap_panel_cf

  poisson_biform(u,v) =  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_cf )dΩ
  poisson_liform(v) = ∫(  (rhs_cf*v)*meas_cf )dΩ
  op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

  # uh = solve(ls,op)

  ## for pvectors, the ghost may not be in the prange of the get_matrix
  ## This causes issues with GridapSolvers Krylov solvers, in the allocation of x
  ## To avoid, allocate x based on the domain of A
  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  uh = FEFunction(U,x)

  e = l2(f_panel_cf-uh,dΩ)

  if return_vtk
    lvl = nref(nc(panel_model))
    cell_geo_map = geo_map_func(Ω_panel)
    panel_cfs = [f_panel_cf,uh,f_panel_cf-uh]
    labels = ["u","uh","eu"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  ### convergence output for DrWatson
  dir_convergence = dir*"/convergence"
  (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

  n = nc(panel_model)
  dxx = dx(panel_model)
  output = @strdict e n dxx p_fe lvl
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("laplace_beltrami_nref$(lvl)_p$p_fe.jld2")), output)

  return e, false,false
end

################################################################################
#### Auto convergence test
################################################################################
function main(distribute,nprocs;octree=false)
  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("Auto conference test: Laplace Beltrami")

  n_ref_lvls = 4
  ps = [1,2,3]
  ls = LUSolver()
  models  = get_refined_models(n_ref_lvls)

  dir = datadir("LaplaceBeltramiConvergence")
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)


  if prod(nprocs) > 1
    i_am_main(ranks) && println("Distributed test")
    if octree
      i_am_main(ranks) && println("Octrees")
      models =  get_octree_refined_models(ranks,n_ref_lvls)
    else
      models,  = get_distributed_refined_models(ranks,nprocs,models)
    end
    ls = CGSolver(JacobiLinearSolver();maxiter=2000,verbose=i_am_main(ranks))
  end

  for (key, val) in analytic_funcs
    i_am_main(ranks) && println("laplace_beltrami_convergence_func_$(key)")
    _dir = dir*"/func_$(key)"
    (i_am_main(ranks) && !isdir(_dir) ) && mkdir(_dir)
    p_convergence_test(ranks,ps,models,laplace_beltrami_solver,_dir,val,ls,true)
  end

  i_am_main(ranks) && println("--DONE--")

end

################################################################################
#### Convergence test with plots
################################################################################

function laplace_beltrami_convergence_test(ranks::AbstractArray,nprocs::Int,
  analytic_funcs,n_ref_lvls=4,ps=[1],ls=LUSolver(),return_vtk=false)

  # serial models
  models  = get_refined_models(n_ref_lvls)

  if prod(nprocs) > 1
    i_am_main(ranks) && println("Distributed test")
    models,  = get_distributed_refined_models(ranks,nprocs,models)
  end

  dir = datadir("LaplaceBeltrami")
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  for (key, val) in analytic_funcs
    simName = "laplace_beltrami_convergence_func_$(key)"
    i_am_main(ranks) && println(simName)

    errors = Vector{Vector{Float64}}(undef,length(ps))
    ns = Vector{Vector{Float64}}(undef,length(ps))
    dxs = Vector{Vector{Float64}}(undef,length(ps))
    slopes = Vector{Float64}(undef,length(ps))

    for (i,p_fe) in enumerate(ps)
      i_am_main(ranks) && println("p_fe = $p_fe")
      errors[i],ns[i],dxs[i],slopes[i] = h_convergence_test(models,laplace_beltrami_solver,p_fe,dir,val,ls,return_vtk)
    end

    i_am_main(ranks) && print_convergence_results(errors,ns,dxs,slopes,ps)

    output = @strdict errors ns dxs slopes ps
    i_am_main(ranks) && safesave(datadir(dir, ("$simName.jld2")), output)

    i_am_main(ranks) && plot_convergence_from_saved(dir,simName)

  end

end



end # module
