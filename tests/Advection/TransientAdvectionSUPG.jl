""" Linear advection equation (material form)
∂ₜu + β⋅ ∇ᵧ(u) = 0
Solve with SUPG as per Brooks & Hughes 1982 paper
Replicate test in Section 2.2, 2.3 of Lauritzen2012 paper
"""

module TransientAdvectionSUPG

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
using CSV
using DataFrames

include("advection_funcs.jl")
include("Lauritzen_functions.jl")
include("../convergence_tools.jl")
include("../output_tools.jl")

using GridapPETSc
function petsc_mumps_setup(ksp)
  pc       = Ref{GridapPETSc.PETSC.PC}()
  mumpsmat = Ref{GridapPETSc.PETSC.Mat}()
  @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPPREONLY)
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCLU)
  @check_error_code GridapPETSc.PETSC.PCFactorSetMatSolverType(pc[],GridapPETSc.PETSC.MATSOLVERMUMPS)
  @check_error_code GridapPETSc.PETSC.PCFactorSetUpMatSolverType(pc[])
  @check_error_code GridapPETSc.PETSC.PCFactorGetMatrix(pc[],mumpsmat)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[],  4, 1)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 28, 2)
  @check_error_code GridapPETSc.PETSC.MatMumpsSetIcntl(mumpsmat[], 29, 1)
  # @check_error_code GridapPETSc.PETSC.MatMumpsSetCntl(mumpsmat[],  1, 0.00001)
  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
end

################################################################################
#### Transient
################################################################################
function transient_advection_supg_solver(
  panel_model::Union{<:DiscreteModel{2,2},<:GridapDistributed.DistributedDiscreteModel{2,2}},
  p_fe::Int,_dir::String,u::Function,v::Function,CFL=0.1,ls=LUSolver(),tF=2*π,return_vtk=false)

  # get the ranks to help with storing/saving solution
  ranks = get_ranks(panel_model)

  lvl = nref(nc(panel_model))
  i_am_main(ranks) && println("nlevl = $lvl")

  dir = _dir*"/sol_p$(p_fe)_nref$lvl"
  (i_am_main(ranks) && !isdir(dir) && return_vtk) && mkdir(dir)

  (i_am_main(ranks) && return_vtk) &&  save_mesh(dir,panel_model)

  ## now enter the solver
  panel_ids = get_panel_ids(panel_model)
  degree = 2*(p_fe + 1)

  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,degree)


  # v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_cf = panelwise_cellfield(u,Ω_panel,panel_ids)

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
  P = TransientTrialFESpace(Q)

  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)

  # supg stabilisation parameter
  _dx = dx(panel_model)
  _dt = _dx*CFL#/p_fe^2
  # dt = floor(_dt,sigdigits=1)
  # dt = _dt

  nsteps = tF/ _dt
  dt = tF/floor(nsteps)

  i_am_main(ranks) && println("nsteps = $nsteps")
  i_am_main(ranks) && println("dt = $dt, other dt = $_dt")

  τ = 0.5*dt

  function get_velocity(t)
    vecX(XYZ) = v(t)(XYZ)
    vX = panel_to_cartesian(tangent_vec(vecX))
    v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
    return v_contr_cf
  end

  β = get_velocity(0.0)
  a_Ω(t, u, v) = ∫( (v * β ⋅ ∇(u))*meas_cf )dΩ
  a_s(t, u, v) = ∫(( (β ⋅ ∇(u)) * (β ⋅ ∇(v)))*meas_cf)dΩ
  m_Ω(t, u, v) = ∫( (∂t(u) * v)*meas_cf)dΩ
  m_s(t, u, v) = ∫( (∂t(u) *(β ⋅ ∇(v)))*meas_cf)dΩ

  m(t, u, v) = m_Ω(t, u, v) + τ * m_s(t, u, v)
  a(t, u, v) = a_Ω(t, u, v) + τ * a_s(t, u, v)

  res(t,u,v) = m(t,u,v) + a(t,u,v)
  jac(t,u,du,v) = a(t,du,v)
  jac_t(t,u,dtu,v) = ∫( v*dtu )dΩ + τ * ∫( dtu *(β ⋅ ∇(v)))dΩ
  opT = TransientFEOperator(res, (jac, jac_t), P, Q)

  # vel = get_velocity(0.0)
  # a_mass_Ω(dtu,v) = ∫( (dtu*v)*meas_cf )dΩ
  # a_mass_s(dtu,v) = ∫( (dtu*(vel⋅∇(v)))*meas_cf )dΩ
  # a_Ω(u,v) = ∫( ((vel⋅∇(u))*v )*meas_cf )dΩ
  # a_s(u,v) =  ∫( ((vel⋅∇(u))*(vel⋅∇(v)) )*meas_cf )dΩ

  # a_mass(t,dtu,v) = a_mass_Ω(dtu,v) + τ*a_mass_s(dtu,v)
  # res(t,u,v) =  a_Ω(u,v) + τ*a_s(u,v)
  # jac(t,u,du,v) = a_Ω(du,v) + τ*a_s(du,v)
  # jac_t(t,u,dtu,v) = a_mass_Ω(dtu,v) + τ*a_mass_s(dtu,v)
  # opT = TransientSemilinearFEOperator(a_mass, res, (jac,jac_t), P, Q, constant_mass=true)

  # a_mass_Ω(dtu,v) = ∫( (dtu*v)*meas_cf )dΩ
  # a_mass_s(t,dtu,v) = ∫( (dtu*(get_velocity(t)⋅∇(v)))*meas_cf )dΩ
  # a_Ω(t,u,v) = ∫( ((get_velocity(t)⋅∇(u))*v )*meas_cf )dΩ
  # a_s(t,u,v) =  ∫( ((get_velocity(t)⋅∇(u))*(get_velocity(t)⋅∇(v)) )*meas_cf )dΩ

  # a_mass(t,dtu,v) = a_mass_Ω(dtu,v) + τ*a_mass_s(t,dtu,v)
  # res(t,u,v) =  a_Ω(t,u,v) + τ*a_s(t,u,v)
  # jac(t,u,du,v) = a_Ω(t,du,v) + τ*a_s(t,du,v)
  # jac_t(t,u,dtu,v) = a_mass_Ω(dtu,v) + τ*a_mass_s(t,dtu,v)
  # opT = TransientSemilinearFEOperator(a_mass, res, (jac,jac_t), P, Q)

  # solve with SSP RK 3
  uh0 = interpolate_everywhere(u_cf, P(0.0))
  # _a(u,v) = ∫( u*v )dΩ
  # _l(v) = ∫( u_cf*v )dΩ
  # op = AffineFEOperator(_a,_l,P,Q)
  # uh0 = solve(LUSolver(),op)
  t0 = 0.0

  nls = GridapSolvers.NonlinearSolvers.NewtonSolver(ls;rtol=1.e-12,verbose=i_am_main(ranks))
  # solver = ThetaMethod(nls,dt,0.5)
  # solver = RungeKutta(nls, ls, dt, :DIRK_CrankNicolson_2_2)
  solver = RungeKutta(nls, ls, dt, :SDIRK_3_2)
  # solver = RungeKutta(nls, ls, dt, :EXRK_SSP_3_3)
  solT = solve(solver, opT, t0, tF, uh0)

  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  cell_geo_map = geo_map_func(Ω_panel)
  labels = ["uh","v", "eu"]

  Ω_error = Triangulation(panel_model)
  dΩ_error = Measure(Ω_error,6*p_fe+1)
  if return_vtk
    panel_cfs = [uh0, covarient_basis_cf⋅ get_velocity(0.0), uh0-uh0]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/solT_0.vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)
    # i_am_main(ranks) && save_cellfields(dir,Ω_panel,t0,[uh0],["uh"])
  end

  ## store errors
  ts = Float64[]
  Es = Float64[]
  Ms = Float64[]
  eu = 0.0
  t = 0.0
  mm = mass_conservation(uh0,meas_cf,dΩ_error)

  push!(ts,0.0)
  push!(Es,0.0)
  push!(Ms,mm)

  counter = 1

  _uh = uh0
  _t = 0.0

  for (t,uh) in solT

    i_am_main(ranks) && println("t = ", t)

    eu = l2((uh-uh0),meas_cf,dΩ_error)
    mm = mass_conservation(uh,meas_cf,dΩ_error)

    push!(ts,t)
    push!(Es,eu)
    push!(Ms,mm)

    _uh = uh
    _t = t

    if return_vtk && (mod(counter,50) == 0)
      panel_cfs = [uh, covarient_basis_cf⋅ get_velocity(t), uh-uh0]
      cellfields = map((x,y) -> x=>y, labels,panel_cfs)

      writevtk(Ω_panel,dir*"/solT_$t.vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)
      # i_am_main(ranks) && save_cellfields(dir,Ω_panel,t,[uh],["uh"])
    end
    counter = counter + 1
  end

  panel_cfs = [_uh, covarient_basis_cf⋅ get_velocity(_t), _uh-uh0]
  cellfields = map((x,y) -> x=>y, labels,panel_cfs)
  writevtk(Ω_panel,_dir*"/solT_final_nref$(lvl)_p$(p_fe).vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)


  push!(ts,_t)
  push!(Es,eu)
  push!(Ms,mm)

  # if return_vtk
  #   if length(ranks) > 1
  #     _make_pvd_distributed(dir,"solT",1)
  #   else
  #     make_pvd(dir,"solT",1)
  #   end
  # end

    ### convergence output for DrWatson
    dir_convergence = _dir*"/convergence"
    (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

    n = nc(panel_model)
    dxx = dx(panel_model)
    output = @strdict n dxx p_fe lvl ts Es Ms
    i_am_main(ranks) && safesave(datadir(dir_convergence, ("transient_advection_dg_nref$(lvl)_p$p_fe.jld2")), output)


  return ts, Es
end

## helper function to return the error for transient solution
function transient_advection_supg_errors(panel_model,args...)
  ts, Es  = transient_advection_supg_solver(panel_model,args...)
  return minimum(Es[end-10:end]),false,false
end


################################################################################
#### Main run for transient solution
################################################################################
function main_transient(distribute;nprocs,options,n_ref_lvls,p_fe,CFL,tF,return_vtk=false)
  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("transient_advection_supg")

  u0(XYZ) = cosine_bell(XYZ)
  u = panel_to_cartesian(u0)
  v = nondivergent_velocity

  parametric_octree_dmodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=n_ref_lvls)
  panel_model = parametric_octree_dmodel.parametric_dmodel
  # panel_model = get_distributed_panel_model(ranks,nprocs,n_ref_lvls)

  dir = datadir("TransientAdvectionSUPG_Laurentz_Octree_$p_fe")
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  GridapPETSc.Init(args=split(options))

  # ls = GMRESSolver(10;Pr=JacobiLinearSolver(),maxiter=2000,verbose=1)
  # ls = LUSolver()
  ls = PETScLinearSolver()

  ts, Es = transient_advection_supg_solver(panel_model,p_fe,dir,u,v,CFL,ls,tF,return_vtk)

  output = @strdict ts Es
  i_am_main(ranks) && safesave(datadir(dir, ("advection_errors.jld2")), output)

  GridapPETSc.Finalize()
  GridapPETSc.gridap_petsc_gc()

  i_am_main(ranks) && println("--DONE--")

end


################################################################################
#### Auto convergence test
################################################################################
function main(distribute,nprocs;octree=false)
  ranks = distribute(LinearIndices((nprocs,)))

  n_ref_lvls = 5
  ps = [3]#[1,2,3]
  ls = LUSolver()
  CFL = 0.5

  v = vt
  u = panel_to_cartesian(u0)
  tF = 2*π

  options = """
  -ksp_type gmres
  -ksp_rtol 1.0e-12
  -ksp_converged_reason
  -ksp_monitor
  """

  dir = foldername("TransientAdvectionSUPG_sdirk",octree,false)
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  models = get_models(ranks,nprocs,n_ref_lvls;threedims=false,octree=octree)

  # GridapPETSc.Init(args=split(options))
  # ls = PETScLinearSolver()

  i_am_main(ranks) && println("transient_advection_supg_convergence")
  p_convergence_test(ranks,ps,models,transient_advection_supg_errors,dir,u,v,CFL,ls,tF,true)

  # GridapPETSc.Finalize()
  # GridapPETSc.gridap_petsc_gc()
end

################################################################################
#### Convergence test with plots
################################################################################
function transient_advection_supg_convergence_test(ranks::AbstractArray,nprocs::Int,
    u::Function,vt::Function,n_ref_lvls=4,ps=[1],CFL=0.1,ls=LUSolver(),tF=2*π,return_vtk=false)

  # serial models
  models  = get_refined_models(n_ref_lvls)

  if prod(nprocs) > 1
    i_am_main(ranks) && println("Distributed test")
    models,  = get_distributed_refined_models(ranks,nprocs,models)
  end

  simName = "transient_advection_supg_convergence"
  i_am_main(ranks) && println(simName)

  dir = datadir("TransientAdvectionSUPG")
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  errors = Vector{Vector{Float64}}(undef,length(ps))
  ns = Vector{Vector{Float64}}(undef,length(ps))
  dxs = Vector{Vector{Float64}}(undef,length(ps))
  slopes = Vector{Float64}(undef,length(ps))


  for (i,p_fe) in enumerate(ps)
    i_am_main(ranks) && println("p_fe = $p_fe")
    errors[i],ns[i],dxs[i],slopes[i] = h_convergence_test(models,transient_advection_supg_errors,p_fe,dir,u,vt,CFL,ls,tF,return_vtk)
  end

  i_am_main(ranks) && print_convergence_results(errors,ns,dxs,slopes,ps)

  output = @strdict errors ns dxs slopes ps
  i_am_main(ranks) && safesave(datadir(dir, ("$simName.jld2")), output)

  i_am_main(ranks) && plot_convergence_from_saved(dir,simName)



end


end # module
