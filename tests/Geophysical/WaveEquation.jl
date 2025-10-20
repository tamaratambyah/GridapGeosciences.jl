"""
solve the linearised wave equation in steady form using manufactured solutions
u + ∇ᵧ(φ) = f₁
φ + ∇ᵧ⋅u = f₁
"""

module WaveEquation

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using GridapGeosciences
using Test

include("../convergence_tools.jl")
include("Williamson2Test.jl")

function wave_solver(panel_model,p_fe::Int,dir::String,h::Function,vX::Function,ls=LUSolver(),return_vtk=false)
  ranks = get_ranks(panel_model)

  lvl = nref(nc(panel_model))
  i_am_main(ranks) && println("nref = $lvl")

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TrialFESpace(V)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
  u_proj_cf = panelwise_cellfield(projection_v(vX),Ω_panel,panel_ids)

  sdiv_cf =  panelwise_cellfield(surfdiv(contra_v(vX)),Ω_panel,panel_ids)

  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  sgrad_cf = panelwise_cellfield(sgrad(h),Ω_panel,panel_ids)
  pinvJ_cf = panelwise_cellfield(forward_pinv_jacobian,Ω_panel,panel_ids)

  # manufacture rhs functions
  rhs_scalar = h_cf + sdiv_cf
  rhs_vector = u_proj_cf + sgrad_cf
  rhs_con_vector = pinvJ_cf ⋅ rhs_vector # exact contravariant component


  metric_cf = CellField(analytic_metric,Ω_panel)
  meas_cf = CellField(sqrtg,Ω_panel)
  grad_meas_cf = CellField(grad_meas,Ω_panel)

  biform1((u,p),(v,q)) = ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
  biform2((u,p),(v,q)) = ∫( (p*q)*meas_cf )dΩ + ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) )  )dΩ

  biformX((u,p),(v,q)) = biform1((u,p),(v,q)) + biform2((u,p),(v,q))
  liformX((v,q)) = ∫( rhs_con_vector⋅(metric_cf⋅v)*meas_cf )dΩ + ∫( (rhs_scalar*q)*meas_cf )dΩ

  op = AffineFEOperator(biformX,liformX,X,Y)
  uh,ph = solve(ls,op)

  uh_proj = covarient_basis_cf ⋅ uh

  e_u = l2( (u_proj_cf - uh_proj)*meas_cf,dΩ) # error in physical velocity u
  e_p = l2((h_cf - ph)*meas_cf,dΩ) # error in depth

  if return_vtk
    lvl = nref(nc(panel_model))
    cell_geo_map = geo_map_func(Ω_panel)
    panel_cfs = [h_cf, u_proj_cf, ph, uh_proj, h_cf-ph, u_proj_cf-uh_proj ]
    labels = ["p","u_proj", "ph", "uh_proj", "ep","eu"]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nref$(lvl)_p$p_fe",cellfields=cellfields,append=false,geo_map=cell_geo_map)

  end

  return e_u,e_p,false
end


function wave_errors(panel_model,p_fe::Int,dir::String,
  h::Function,vX::Function,f::Function,η::Function,b::Function,ls=LUSolver(),CFL=0.1,return_vtk=false)
  wave_solver(panel_model,p_fe,dir,h,vX,ls,return_vtk)
end


################################################################################
#### Auto convergence test
################################################################################
function main(distribute,nprocs;octree=false)
  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("Auto conference test: Wave Equation")

  n_ref_lvls = 4
  ps = [1,2,3]
  ζs = [0.0]
  ls = LUSolver()
  models  = get_refined_models(n_ref_lvls)

  dir = datadir("WaveConvergence")
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  if prod(nprocs) > 1
    i_am_main(ranks) && println("Distributed test")
    if octree
      i_am_main(ranks) && println("Octrees")
      models =  get_octree_refined_models(ranks,n_ref_lvls)
    else
      models,  = get_distributed_refined_models(ranks,nprocs,models)
    end
    # ls = CGSolver(JacobiLinearSolver();maxiter=2000,verbose=i_am_main(ranks))
  end



  for (i,ζ) in enumerate(ζs)
    _dir = dir*"/func_z$i"
    (i_am_main(ranks) && !isdir(_dir) ) && mkdir(_dir)

    h = panel_to_cartesian(h₀(ζ))
    vX = panel_to_cartesian(tangent_vec(u₀(ζ)))

    i_am_main(ranks) && println("wave_equation_convergence_func_z$i")
    p_convergence_test(ranks,ps,models,wave_solver,_dir,h,vX,ls,true)
  end

  i_am_main(ranks) && println("--DONE--")
end

################################################################################
#### Convergence test with plots
################################################################################

function wave_convergence_test(ranks::AbstractArray,nprocs::Int,
  ζs=[0.0],n_ref_lvls=4,ps=[1],ls=LUSolver(),CFL=0.1,return_vtk=false)

  williamson2_convergence_test(ranks,nprocs,wave_errors,ζs,n_ref_lvls,ps,ls,CFL,return_vtk)
end

end # module
