""" Linear advection equation (material form)
∂ₜu + β⋅ ∇ᵧ(u) = 0
Solve with SUPG as per Brooks & Hughes 1982 paper
Replicate test in Section 2.2, 2.3 of Lauritzen2012 paper
"""

module TransientAdvectionSUPG_Lauentz

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

include("Lauritzen_functions.jl")
include("../convergence_tools.jl")
include("../output_tools.jl")

using GridapPETSc

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
  uh0 = interpolate_everywhere(u_cf, P(0.0))
  t0 = 0.0


  nls = GridapSolvers.NonlinearSolvers.NewtonSolver(ls;rtol=1.e-14,verbose=i_am_main(ranks))
  solver = RungeKutta(nls, ls, dt, :DIRK_CrankNicolson_2_2)
  # solver = RungeKutta(ls, ls, dt, :EXRK_SSP_3_3)
  solT = solve(solver, opT, t0, tF, uh0)

  covariant_basis_cf = panelwise_cellfield(covariant_basis,Ω_panel,panel_ids)

  cell_geo_map = geo_map_func(Ω_panel)
  labels = ["uh","v", "eu"]

  Ω_error = Triangulation(panel_model)
  dΩ_error = Measure(Ω_error,6*p_fe+1)
  if return_vtk
    panel_cfs = [uh0, covariant_basis_cf⋅ get_velocity(0.0), uh0-uh0]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk_with_cell_geomap(cell_geo_map,Ω_panel,dir*"/solT_0.vtu", cellfields=cellfields,append=false)
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
      panel_cfs = [uh, covariant_basis_cf⋅ get_velocity(t), uh-uh0]
      cellfields = map((x,y) -> x=>y, labels,panel_cfs)

      writevtk_with_cell_geomap(cell_geo_map,Ω_panel,dir*"/solT_$t.vtu", cellfields=cellfields,append=false)
      # i_am_main(ranks) && save_cellfields(dir,Ω_panel,t,[uh],["uh"])
    end
    counter = counter + 1
  end

  panel_cfs = [_uh, covariant_basis_cf⋅ get_velocity(_t), _uh-uh0]
  cellfields = map((x,y) -> x=>y, labels,panel_cfs)
  writevtk_with_cell_geomap(cell_geo_map,Ω_panel,_dir*"/solT_final_nref$(lvl)_p$(p_fe).vtu", cellfields=cellfields,append=false)


  push!(ts,_t)
  push!(Es,eu)
  push!(Ms,mm)


  ### convergence output for DrWatson
  dir_convergence = _dir*"/convergence"
  (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

  n = nc(panel_model)
  dxx = dx(panel_model)
  output = @strdict n dxx p_fe lvl ts Es Ms
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("transient_advection_dg_nref$(lvl)_p$p_fe.jld2")), output)


  return ts, Es
end




################################################################################
#### Main run for transient solution
################################################################################
function main_transient(distribute;nprocs,options,n_ref_lvls,p_fe,CFL,tF,return_vtk=false)
  ranks = distribute(LinearIndices((nprocs,)))

  radius = 1.0

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("transient_advection_supg")

  u0(XYZ) = cosine_bell(XYZ)
  u = panel_to_cartesian(u0)
  v = nondivergent_velocity

  parametric_octree_dmodel = ParametricOctreeDistributedDiscreteModel(ranks, radius; num_initial_uniform_refinements=n_ref_lvls)
  panel_model = parametric_octree_dmodel.parametric_dmodel

  dir = datadir("TransientAdvectionSUPG_Laurentz_Octree_$(p_fe)_CN")
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







end # module
