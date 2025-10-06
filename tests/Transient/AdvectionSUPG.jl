module AdvectionSUPG

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

include("Lauritzen_functions.jl")
include("../helpers.jl")
################################################################################
#### Transient
################################################################################
function transient_advection_supg(panel_model,p_fe::Int,u::Function,v::Function,CFL=0.1,ls=LUSolver(),tF=2*π,return_vtk=false)
  lvl = nref(nc(panel_model))
  println("nlevl = $lvl")

  dir = datadir("Transient_advection_nref$lvl")
  (return_vtk && !isdir(dir)) && mkdir(dir)


  panel_ids = get_panel_ids(panel_model)
  degree = 2*(p_fe + 1)

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)


  # v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_cf = panelwise_cellfield(u,Ω_panel,panel_ids)

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  P = TrialFESpace(Q)

  # hard code RT space as order 1 -- for velocity
  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
  U = TrialFESpace(V)

  meas_cf = CellField(sqrtg,Ω_panel)

  # supg stabilisation parameter
  _dx = dx(nc(panel_model))
  _dt = _dx*CFL/p_fe
  dt = floor(_dt,sigdigits=1)
  τ = 0.5*dt

  function get_velocity(t)
    vecX(XYZ) = v(t)(XYZ)
    vX = panel_to_cartesian(tangent_vec(vecX))
    v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
    return v_contr_cf
    # interpolate(v_contr_cf,U)
  end


  # a_mass_Ω(dtu,v) = ∫( (dtu*v)*meas_cf )dΩ
  # a_mass_s(dtu,v) = ∫( (dtu*(vel⋅∇(v)))*meas_cf )dΩ
  # a_Ω(u,v) = ∫( ((vel⋅∇(u))*v )*meas_cf )dΩ
  # a_s(u,v) =  ∫( ((vel⋅∇(u))*(vel⋅∇(v)) )*meas_cf )dΩ

  # a_mass(t,dtu,v) = a_mass_Ω(dtu,v) + τ*a_mass_s(dtu,v)
  # res(t,u,v) =  a_Ω(u,v) + τ*a_s(u,v)
  # jac(t,u,du,v) = a_Ω(du,v) + τ*a_s(du,v)
  # jac_t(t,u,dtu,v) = a_mass_Ω(dtu,v) + τ*a_mass_s(dtu,v)
  # opT = TransientSemilinearFEOperator(a_mass, res, (jac,jac_t), P, Q, constant_mass=true)

  a_mass_Ω(dtu,v) = ∫( (dtu*v)*meas_cf )dΩ
  a_mass_s(t,dtu,v) = ∫( (dtu*(get_velocity(t)⋅∇(v)))*meas_cf )dΩ
  a_Ω(t,u,v) = ∫( ((get_velocity(t)⋅∇(u))*v )*meas_cf )dΩ
  a_s(t,u,v) =  ∫( ((get_velocity(t)⋅∇(u))*(get_velocity(t)⋅∇(v)) )*meas_cf )dΩ

  a_mass(t,dtu,v) = a_mass_Ω(dtu,v) + τ*a_mass_s(t,dtu,v)
  res(t,u,v) =  a_Ω(t,u,v) + τ*a_s(t,u,v)
  jac(t,u,du,v) = a_Ω(t,du,v) + τ*a_s(t,du,v)
  jac_t(t,u,dtu,v) = a_mass_Ω(dtu,v) + τ*a_mass_s(t,dtu,v)
  opT = TransientSemilinearFEOperator(a_mass, res, (jac,jac_t), P, Q)

  # solve with SSP RK 3
  uh0 = interpolate(u_cf, P)
  t0 = 0.0

  solver = RungeKutta(ls, ls, dt, :EXRK_SSP_3_3)
  solT = solve(solver, opT, t0, tF, uh0)

  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  cell_geo_map = geo_map_func(panel_ids)

  if return_vtk
    # writevtk(Ω_panel,dir*"/solT_0.vtu", cellfields=["uh"=>uh0,"v"=>covarient_basis_cf⋅ vel],append=false,geo_map=cell_geo_map)
    writevtk(Ω_panel,dir*"/solT_0.vtu", cellfields=["uh"=>uh0,"v"=>covarient_basis_cf⋅ get_velocity(0.0)],append=false,geo_map=cell_geo_map)
  end

  ## store errors
  ts = Float64[]
  Es = Float64[]

  push!(ts,0.0)
  push!(Es,0.0)

  for (t,uh) in solT

    println(t)

    eu = l2((uh-uh0)*meas_cf,dΩ)

    push!(ts,t)
    push!(Es,eu)
    if return_vtk
      # writevtk(Ω_panel,dir*"/solT_$t.vtu", cellfields=["uh"=>uh],append=false,geo_map=cell_geo_map)
      writevtk(Ω_panel,dir*"/solT_$t.vtu", cellfields=["uh"=>uh,"v"=>covarient_basis_cf⋅ get_velocity(t)],append=false,geo_map=cell_geo_map)
    end
  end

  output = @strdict ts Es
  safesave(datadir(dir, ("advection_errors.jld2")), output)

  return Es[end],false,false
end

################################################################################
#### Main run for transient solution
################################################################################
function main(distribute;nprocs,options,n_ref_lvls,p_fe,CFL,tF,return_vtk=false)
  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("advection_supg_convergence_func")

  u0(XYZ) = cosine_bell(XYZ)
  u = panel_to_cartesian(u0)
  v = nondivergent_velocity

  GridapPETSc.Init(args=split(options))

  # ls = GMRESSolver(10;Pr=JacobiLinearSolver(),maxiter=2000,verbose=1)
  # ls = LUSolver()
  ls = PETScLinearSolver()

  panel_model = get_distributed_panel_model(ranks,nprocs,n_ref_lvls)

  transient_advection_supg(panel_model,p_fe,u,v,CFL,ls,tF,return_vtk)

  GridapPETSc.Finalize()
  GridapPETSc.gridap_petsc_gc()

  i_am_main(ranks) && println("--DONE--")

end


################################################################################
#### Convergence test with plots
################################################################################
# function transient_advection_supg_convergence_test(ranks::AbstractArray,nprocs::Int,
#   dir,u::Function,vX::Function,n_ref_lvls=4,ps=[1],CFL=0.1,ls=LUSolver(),tF=2*π,return_vtk=false)

#   # serial models
#   models  = get_refined_models(n_ref_lvls)

#   if prod(nprocs) > 1
#     i_am_main(ranks) && println("Distributed test")
#     models,  = get_distributed_refined_models(ranks,nprocs,models)
#   end

#   simName = "transient_advection_supg_convergence_func"
#   i_am_main(ranks) && println(simName)

#   errors = Vector{Vector{Float64}}(undef,length(ps))
#   ns = Vector{Vector{Float64}}(undef,length(ps))
#   dxs = Vector{Vector{Float64}}(undef,length(ps))
#   slopes = Vector{Float64}(undef,length(ps))


#   for (i,p_fe) in enumerate(ps)
#     i_am_main(ranks) && println("p_fe = $p_fe")
#     errors[i],ns[i],dxs[i],slopes[i] = h_convergence_test(models,advection_supg_solver,p_fe,u,vX,CFL,ls,tF,return_vtk)
#   end

#   i_am_main(ranks) && print_convergence_results(errors,ns,dxs,slopes,ps)

#   output = @strdict errors ns dxs slopes ps
#   i_am_main(ranks) && safesave(datadir(dir, ("$simName.jld2")), output)

#   i_am_main(ranks) && plot_convergence_from_saved(dir,simName)

#   if return_vtk
#     if nprocs > 1
#       _make_pvd_distributed(dir,"solT",1)
#     else
#       make_pvd(dir,"solT",1)
#     end
#   end

# end


end # module
