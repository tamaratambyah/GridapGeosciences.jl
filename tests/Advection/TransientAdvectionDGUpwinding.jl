""" Linear advection equation (flux form)
∂ₜu + ∇ᵧ⋅(βu) = 0
Solve with dG upwinding as per Brezzi 2004 paper
Replicate test in Section 5.4 of Rognes2013 paper
"""

module TransientAdvectionDGUpwinding

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using MPIPreferences
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra
using GridapGeosciences
using GridapPETSc
using Test

include("advection_funcs.jl")
include("../convergence_tools.jl")

function my_mean( Bu_n::SkeletonPair)
  plus  = ( Bu_n.plus)
  minus = ( Bu_n.minus)
  0.5*( plus - minus  )
end

function transient_advection_dg_solver(
  panel_model::Union{<:DiscreteModel{2,2},<:GridapDistributed.DistributedDiscreteModel{2,2}},
  p_fe::Int,_dir::String,u::Function,vX::Function,CFL=0.1,ls=LUSolver(),tF=2*π,return_vtk=false)

  # get the ranks to help with storing/saving solution
  ranks = get_ranks(panel_model)

  lvl = nref(nc(panel_model))
  i_am_main(ranks) && println("nlevl = $lvl")

  dir = _dir*"/sol_p$(p_fe)_nref$lvl"
  (i_am_main(ranks) && !isdir(dir) && return_vtk) && mkdir(dir)

  panel_ids = get_panel_ids(panel_model)
  degree = 2*(p_fe + 1)

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)

  Λ = SkeletonTriangulation(panel_model)
  if typeof(Ω_panel) <: GridapDistributed.DistributedTriangulation
    i_am_main(ranks) && println("ghost skel mesh")
    Λ = SkeletonTriangulation(with_ghost,panel_model)
  end
  dΛ = Measure(Λ,degree)
  n_Λ = get_normal_vector(Λ)

  v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_cf = panelwise_cellfield(u,Ω_panel,panel_ids)

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TransientTrialFESpace(Q)

  # hard code RT space as order 1 -- for velocity
  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
  U = TrialFESpace(V)

  # initial conditions
  # vel = interpolate(v_contr_cf,U)
  vel = v_contr_cf
  uh0 = interpolate(u_cf, P(0.0))
  # _a(u,v) = ∫( u*v )dΩ
  # _l(v) = ∫( u_cf*v )dΩ
  # op = AffineFEOperator(_a,_l,P(0.0),Q)
  # uh0 = solve(LUSolver(),op)

  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  meas_cf_skel = panelwise_cellfield(sqrtg,Λ)

  ## weak form
  a_mass(t,dtu,v) = ∫( (dtu*v)*meas_cf )dΩ

  a_Ω(u,v) =   ∫( -(u*(∇(v)⋅vel) )*meas_cf )dΩ
  a_s1(u,v) = ∫( my_mean((vel*u)⋅n_Λ)*jump(v)*meas_cf_skel.plus   )dΛ

  upwind = abs((vel⋅ n_Λ).plus)/2
  a_s2(u,v) = ∫(  upwind*jump(u)*jump(v)*meas_cf_skel.plus   )dΛ

  res(t,u,v) =  a_Ω(u,v) + a_s1(u,v) + a_s2(u,v)
  jac(t,u,du,v) = a_Ω(du,v) + a_s1(du,v) + a_s2(du,v)
  jac_t(t,u,dtu,v) = ∫( (dtu*v)*meas_cf )dΩ
  opT = TransientSemilinearFEOperator(a_mass, res, (jac,jac_t), P, Q, constant_mass=true)

  # solve with SSP RK 3
  t0 = 0.0
  _dt = dx(panel_model)*CFL/p_fe
  dt = floor(_dt,sigdigits=1)


  solver = RungeKutta(ls, ls, dt, :EXRK_SSP_3_3)
  solT = solve(solver, opT, t0, tF, uh0)

  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  cell_geo_map = geo_map_func(Ω_panel)
  if return_vtk
    writevtk(Ω_panel,dir*"/solT_0.vtu", cellfields=["uh"=>uh0,"v"=>covarient_basis_cf⋅ vel],append=false,geo_map=cell_geo_map)
  end

  ## store errors
  ts = Float64[]
  Es = Float64[]

  counter = 1
  eu = 0.0
  t = 0.0

  push!(ts,t)
  push!(Es,eu)

  Ω_error = Triangulation(panel_model)
  dΩ_error = Measure(Ω_error,6*p_fe+1)
  for (t,uh) in solT

    i_am_main(ranks) && println(t)

    eu = l2((uh-uh0),meas_cf,dΩ_error)

    push!(ts,t)
    push!(Es,eu)

    if return_vtk  && (mod(counter,10) == 0)
      writevtk(Ω_panel,dir*"/solT_$t.vtu", cellfields=["uh"=>uh],append=false,geo_map=cell_geo_map)
    end
    counter = counter + 1
  end

  push!(ts,t)
  push!(Es,eu)

  if return_vtk
    if length(ranks) > 1
      _make_pvd_distributed(dir,"solT",1)
    else
      make_pvd(dir,"solT",1)
    end
  end



  return ts, Es

end


function transient_advection_dg_errors(panel_model,args...)
  ts, Es  = transient_advection_dg_solver(panel_model,args...)
  return minimum(Es[end-10:end]),false,false
end

################################################################################
#### Auto convergence test
################################################################################
function main(distribute,nprocs;octree=false)
  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("Auto conference test: Transient AdvectionDGUpwinding")

  n_ref_lvls = 4
  ps = [1,2]#[1,2,3]
  ls = LUSolver()
  CFL = 0.1

  v = panel_to_cartesian(tangent_vec(vecX))
  u = panel_to_cartesian(u0)
  tF = 2*π

  models  = get_refined_models(n_ref_lvls)

  dir = datadir("TransientAdvectionDGUpwinding")
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


  i_am_main(ranks) && println("transient_advection_dg_convergence")
  p_convergence_test(ranks,ps,models,transient_advection_dg_errors,dir,u,v,CFL,ls,tF,true)

  i_am_main(ranks) && println("--DONE--")
end

################################################################################
#### Convergence test with plots
################################################################################
function transient_advection_dg_convergence_test(ranks::AbstractArray,nprocs::Int,
  u::Function,vX::Function,n_ref_lvls=4,ps=[1],CFL=0.1,ls=LUSolver(),tF=2*π,return_vtk=false)

# serial models
models  = get_refined_models(n_ref_lvls)

if prod(nprocs) > 1
  i_am_main(ranks) && println("Distributed test")
  models,  = get_distributed_refined_models(ranks,nprocs,models)
end

simName = "transient_advection_dg_convergence"
i_am_main(ranks) && println(simName)

dir = datadir("TransientAdvectionDG")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

errors = Vector{Vector{Float64}}(undef,length(ps))
ns = Vector{Vector{Float64}}(undef,length(ps))
dxs = Vector{Vector{Float64}}(undef,length(ps))
slopes = Vector{Float64}(undef,length(ps))


for (i,p_fe) in enumerate(ps)
  i_am_main(ranks) && println("p_fe = $p_fe")
  errors[i],ns[i],dxs[i],slopes[i] = h_convergence_test(models,transient_advection_dg_errors,p_fe,dir,u,vX,CFL,ls,tF,return_vtk)
end

i_am_main(ranks) && print_convergence_results(errors,ns,dxs,slopes,ps)

output = @strdict errors ns dxs slopes ps
i_am_main(ranks) && safesave(datadir(dir, ("$simName.jld2")), output)

i_am_main(ranks) && plot_convergence_from_saved(dir,simName)



end


end ## module
