""" Linear advection equation (material form)
∂ₜu + β⋅ ∇ᵧ(u) = 0
Solve with SUPG as per Brooks & Hughes 1982 paper
Replicate test in Section 5.4 of Rognes2013 paper
"""

module AdvectionSUPG

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using GridapGeosciences
using Test


include("advection_funcs.jl")
include("../convergence_tools.jl")

################################################################################
#### Steady with manufactured solutions
################################################################################
function advection_supg_solver(panel_model,p_fe::Int,dir,u::Function,vX::Function,CFL=0.1,ls=LUSolver(),return_vtk=false)
  lvl = nref(nc(panel_model))
  println("nref = $lvl")

  panel_ids = get_panel_ids(panel_model)
  degree = 2*(p_fe + 1)

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)

  _rhs(p) = αβ -> u(p)(αβ) + vX(p)(αβ)⋅sgrad(u,p)(αβ)

  v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_cf = panelwise_cellfield(u,Ω_panel,panel_ids)
  rhs_cf = panelwise_cellfield(_rhs,Ω_panel,panel_ids)

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  P = TrialFESpace(Q)

  # hard code RT space as order 1 -- for velocity
  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
  U = TrialFESpace(V)

  # vel = interpolate(v_contr_cf,U)
  _a(u,v) = ∫( u⋅v )dΩ
  _l(v) = ∫( v_contr_cf ⋅ v )dΩ
  op = AffineFEOperator(_a,_l,U,V)
  vel = solve(LUSolver(),op)


  meas_cf = CellField(sqrtg,Ω_panel)

  # supg stabilisation parameter
  _dx = dx(nc(panel_model))
  _dt = _dx*CFL/p_fe
  dt = floor(_dt,sigdigits=1)
  τ = 0.5*dt

  a_Ω(u,v) = ∫( (u*v)*meas_cf )dΩ + ∫( ((vel⋅∇(u))*v )*meas_cf )dΩ
  a_s(u,v) =  ∫( (u*(vel⋅∇(v)) )*meas_cf )dΩ + ∫( ((vel⋅∇(u))*(vel⋅∇(v)) )*meas_cf )dΩ

  l_Ω(v) = ∫( rhs_cf*v*meas_cf )dΩ
  l_s(v) = ∫( rhs_cf*(vel⋅∇(v))*meas_cf )dΩ

  biform_advection(u,v) = a_Ω(u,v) + τ*a_s(u,v)
  liform_advection(v) = l_Ω(v) + τ*l_s(v)

  op = AffineFEOperator(biform_advection,liform_advection,P,Q)

  # uh = solve(ls,op)
  A = get_matrix(op)
  b = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b)
  uh = FEFunction(P,x)

  eu = l2((uh-u_cf)*meas_cf,dΩ)

  if return_vtk
    cell_geo_map = geo_map_func(Ω_panel)
    labels = ["uh","u","eu"]
    panel_cfs = [uh,u_cf,uh-u_cf]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe", cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  return eu,false,false
end

################################################################################
#### Auto convergence test
################################################################################
function main(distribute,nprocs)
  ranks = distribute(LinearIndices((nprocs,)))

  n_ref_lvls = 4
  ps = [1,2,3]
  ls = LUSolver()
  CFL = 0.1

  vX = panel_to_cartesian(tangent_vec(vecX))
  u = panel_to_cartesian(u0)

  models  = get_refined_models(n_ref_lvls)

  if prod(nprocs) > 1
    i_am_main(ranks) && println("Distributed test")
    models,  = get_distributed_refined_models(ranks,nprocs,models)
    # ls = CGSolver(JacobiLinearSolver();maxiter=2000,verbose=i_am_main(ranks))
  end

  i_am_main(ranks) && println("advection_supg_convergence_func")
  p_convergence_test(ranks,ps,models,advection_supg_solver,"",u,vX,CFL,ls)

end


################################################################################
#### Convergence test with plots
################################################################################
function advection_supg_convergence_test(ranks::AbstractArray,nprocs::Int,
  u::Function,vX::Function,n_ref_lvls=4,ps=[1],CFL=0.1,ls=LUSolver(),return_vtk=false)

  # serial models
  models  = get_refined_models(n_ref_lvls)

  if prod(nprocs) > 1
    i_am_main(ranks) && println("Distributed test")
    models,  = get_distributed_refined_models(ranks,nprocs,models)
  end

  simName = "advection_supg_convergence_func"
  i_am_main(ranks) && println(simName)

  dir = datadir("AdvectionSUPG")
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  errors = Vector{Vector{Float64}}(undef,length(ps))
  ns = Vector{Vector{Float64}}(undef,length(ps))
  dxs = Vector{Vector{Float64}}(undef,length(ps))
  slopes = Vector{Float64}(undef,length(ps))

  for (i,p_fe) in enumerate(ps)
    i_am_main(ranks) && println("p_fe = $p_fe")
    errors[i],ns[i],dxs[i],slopes[i] = h_convergence_test(models,advection_supg_solver,p_fe,dir,u,vX,CFL,ls,return_vtk)
  end

  i_am_main(ranks) && print_convergence_results(errors,ns,dxs,slopes,ps)

  output = @strdict errors ns dxs slopes ps
  i_am_main(ranks) && safesave(datadir(dir, ("$simName.jld2")), output)

  i_am_main(ranks) && plot_convergence_from_saved(dir,simName)


end


end
