"""
solve the linearised wave equation in
∂ₜu + ∇ᵧ(φ) = 0
∂ₜφ + ∇ᵧ⋅u = 0
"""


module TransientWaveEquation

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


include("../convergence_tools.jl")

## initial conditions
vecX(XYZ) = zero(XYZ)
depth(XYZ) = 1.0 + 0.01*exp(-5*((1-XYZ[1])^2+(0-XYZ[2])^2+(0-XYZ[3])^2))


function transient_wave_solver(panel_model,p_fe::Int,_dir::String,
    h::Function,vX::Function,CFL=0.1,ls=LUSolver(),tF=2*π,return_vtk=false)

  # get the ranks to help with storing/saving solution
  ranks = get_ranks(panel_model)

  lvl = nref(nc(panel_model))
  i_am_main(ranks) && println("nlevl = $lvl")

  dir = _dir*"/sol_p$(p_fe)_nref$lvl"
  (i_am_main(ranks) && !isdir(dir) && return_vtk) && mkdir(dir)


  ## finite element solver
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TransientTrialFESpace(Q)

  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TransientTrialFESpace(V)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  ## initial conditions
  h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
  vec_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  xh0 = interpolate([vec_contra_cf,h_cf],X)
  # _a((u,p),(v,q)) = ∫( u⋅v + p*q )dΩ
  # _l(v) = ∫( vec_contra_cf⋅v + h_cf*q )dΩ
  # op = AffineFEOperator(_a,_l,X,Y)
  # xh0 = solve(LUSolver(),op)

  ## transient weak form
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)

  mass(t, (dtu,dtp), (v,q)) = ∫( (v⋅ (metric_cf⋅ dtu) )*meas_cf )dΩ + ∫( (q*dtp)*meas_cf )dΩ
  res(t,(u,p),(v,q)) =  ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) )  )dΩ - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
  opT = TransientSemilinearFEOperator(mass, res, X, Y, constant_mass=true)


  # transient parameters
  t0 = 0.0
  _dt = dx(nc(panel_model))*CFL/p_fe
  dt = floor(_dt, sigdigits=1)

  # solve with SSP RK 3
  solver = RungeKutta(ls, ls, dt, :EXRK_SSP_3_3)
  solT = solve(solver, opT, t0, tF, xh0)

  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  cell_geo_map = geo_map_func(Ω_panel)

  ## iterate solution
  labels = ["uh","ph"]
  if return_vtk
    panel_cfs = [covarient_basis_cf⋅xh0[1], xh0[2]]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/solT_0" * ".vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  ## casimirs to store
  ts = Float64[]
  ms = Float64[]
  Es = Float64[]
  s_divus = Float64[]
  divus = Float64[]


  _m0 = sum( ∫(  meas_cf*xh0[2] )dΩ)
  _E0 = sum( ∫( 0.5*( xh0[1]⋅(metric_cf⋅ xh0[1]) + xh0[2]*xh0[2])*meas_cf )dΩ  )
  _s_divu = sum(∫(  divergence(meas_cf*xh0[1]) )dΩ)
  _divu = sum(∫(  divergence(xh0[1]) )dΩ)

  push!(ts,0)
  push!(ms,_m0)
  push!(Es,_E0)
  push!(s_divus,_s_divu)
  push!(divus,_divu)

  counter = 1
  for (t, xh) in solT
    uh,ph = xh
    i_am_main(ranks) && println("t = ", t)

    _m = sum( ∫(  meas_cf*ph )dΩ)
    _E = sum( ∫( 0.5*( uh⋅(metric_cf⋅ uh) + ph*ph)*meas_cf )dΩ  )
    _s_divu = sum(∫(   divergence(meas_cf*uh) )dΩ)
    _divu = sum(∫(  divergence(uh)  )dΩ)


    push!(ts,t)
    push!(ms,_m)
    push!(Es,_E)
    push!(s_divus,_s_divu)
    push!(divus,_divu)

    if return_vtk  && (mod(counter,10) == 0)
      panel_cfs = [covarient_basis_cf⋅uh, ph]
      cellfields = map((x,y) -> x=>y, labels,panel_cfs)
      writevtk(Ω_panel,dir*"/solT_$t" * ".vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)
    end
    counter = counter + 1
  end


  if return_vtk
    if length(ranks) > 1
      _make_pvd_distributed(dir,"solT",1)
    else
      make_pvd(dir,"solT",1)
    end
  end
  return ts, Es, ms, s_divus, divus
end


################################################################################
#### Main run for transient solution
################################################################################
function main_transient(distribute,nprocs;options="",n_ref_lvls=4,p_fe=1,CFL=0.1,tF=2*π,return_vtk=false)
  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("transient_wave_equation")

  h = panel_to_cartesian(depth)
  vX = panel_to_cartesian(tangent_vec(vecX))
  ls = LUSolver()

  models  = get_refined_models(n_ref_lvls)

  if prod(nprocs) > 1
    i_am_main(ranks) && println("Distributed test")
    models,  = get_distributed_refined_models(ranks,nprocs,models)
  end

  panel_model = models[1]

  dir = datadir("Transient_wave_equation")
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  GridapPETSc.Init(args=split(options))

  ls = CGSolver(JacobiLinearSolver();maxiter=2000,verbose=i_am_main(ranks))
  # ls = LUSolver()
  # ls = PETScLinearSolver()

  ts, Es, ms, s_divus, divus = transient_wave_solver(panel_model,p_fe,dir,h,vX,CFL,ls,tF,return_vtk)

  output = @strdict ts Es ms s_divus divus
  i_am_main(ranks) && safesave(datadir(dir, ("wave_errors.jld2")), output)

  GridapPETSc.Finalize()
  GridapPETSc.gridap_petsc_gc()

  i_am_main(ranks) && println("--DONE--")
  @test true
end



end ##module
