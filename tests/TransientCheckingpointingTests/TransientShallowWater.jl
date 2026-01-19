"""
solve the non-linear shallow water equations
∂ₜu + q F^† + ∇ᵧ(Φ) = 0
∂ₜφ + ∇ᵧ⋅F = 0
F = φu
Φ = 0.5(u⋅u) + gᵣφ
q = 1/φ( ∇ᵧ^†⋅u  + f )
"""

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using MPI
using PartitionedArrays
using MPIPreferences
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra
using GridapGeosciences
using GridapPETSc
using Test

include("../convergence_tools.jl")
include("helpers.jl")
include("../Geophysical/Williamson2Test.jl")
# include("../Geophysical/Williamson5Test.jl")



function transient_shallow_water_solver(panel_model::Union{<:DiscreteModel{2,2},<:GridapDistributed.DistributedDiscreteModel{2,2}},
  p_fe::Int,dir::String,h::Function,vX::Function,f::Function,b::Function,
  CFL=0.1,lss=(LUSolver(),LUSolver()),restart=false)

  ls_ode, ls_diag = lss

  # get the ranks to help with storing/saving solution
  ranks = get_ranks(panel_model)

  sim_dir = dir*"/sim_data"
  (i_am_main(ranks) && !isdir(sim_dir) ) && mkdir(sim_dir)

  final_dir = dir*"/final_solution"
  (i_am_main(ranks) && !isdir(final_dir) ) && mkdir(final_dir)

  initial_dir = dir*"/initial_solution"
  (i_am_main(ranks) && !isdir(initial_dir) ) && mkdir(initial_dir)

  prog_dir = sim_dir*"/prognostics"
  (i_am_main(ranks) && !isdir(prog_dir) ) && mkdir(prog_dir)

  diag_dir = sim_dir*"/diagnostics"
  (i_am_main(ranks) && !isdir(diag_dir) ) && mkdir(diag_dir)


  ## finite element solver
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  R = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1)
  H = TransientTrialFESpace(R)

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TransientTrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TransientTrialFESpace(V)

  X_prog = TransientMultiFieldFESpace([U,P]) # u, p
  Y_prog = MultiFieldFESpace([V,Q]) # u, p

  X_diag = TransientMultiFieldFESpace([H,U,P]) # q, F, Φ
  Y_diag = MultiFieldFESpace([R,V,Q]) # q, F, Φ

  ## initial conditions
  function initial_condition()
    i_am_main(ranks) && println("initial condition")

    u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
    u_contra_h = interpolate(u_contra_cf,U)

    h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
    b_cf = panelwise_cellfield(b,Ω_panel,panel_ids)
    h_h = interpolate(h_cf-b_cf,P)

    xh0 = interpolate_everywhere([u_contra_h,h_h],X_prog(0.0))
    t = 0.0
    psave(prog_dir*"/solT_$(t)",xh0)
    psave(initial_dir*"/solT_$(t)",xh0)
    return t,xh0
  end

  simName = "solT"
  t0,xh0 = (restart) ? load_last(ranks,X_prog(0.0),prog_dir,simName) : initial_condition()

  ## transient weak form
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  cor_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  gravity = _g
  b_cf = panelwise_cellfield(b,Ω_panel,panel_ids)

  #### DIAGNOSTIC VARIABLES
  # vorticity
  perp_matrix_cf = panelwise_cellfield(perp_matrix,Ω_panel,panel_ids)
  resq(((u,p),(q,F,Φ)),(w,v,ψ)) = ∫( q*p*w*meas_cf  )dΩ - ∫( cor_cf*w*meas_cf  )dΩ - ∫( (perp_matrix_cf⋅u)⋅∇(w)  )dΩ

  # mass flux
  resF(((u,p),(q,F,Φ)),(w,v,ψ)) = ∫( (F⋅ (metric_cf⋅v))*meas_cf )dΩ - ∫( p*(u⋅(metric_cf⋅v))*meas_cf   )dΩ

  # Bernoulli potential
  resΦ(((u,p),(q,F,Φ)),(w,v,ψ)) = ∫( Φ*ψ*meas_cf  )dΩ - ∫( gravity*(p+b_cf)*ψ*meas_cf  )dΩ - ∫( 0.5*( u ⋅(metric_cf⋅u) )ψ*meas_cf  )dΩ

  res_y(t,((u,p),(q,F,Φ)),(w,v,ψ)) = resq(((u,p),(q,F,Φ)),(w,v,ψ)) + resF(((u,p),(q,F,Φ)),(w,v,ψ)) + resΦ(((u,p),(q,F,Φ)),(w,v,ψ))
  jac_y(t,((u,p),(q,F,Φ)),(dq,dF,dΦ),(w,v,ψ)) = ∫( dq*p*w*meas_cf  )dΩ + ∫( (dF⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( dΦ*ψ*meas_cf  )dΩ

  _res_y((q,F,Φ),(w,v,ψ))  = res_y(t0,(xh0,(q,F,Φ)),(w,v,ψ))
  _jac_y((q,F,Φ),(dq,dF,dΦ),(w,v,ψ)) = jac_y(t0,(xh0,(q,F,Φ)),(dq,dF,dΦ),(w,v,ψ))
  _opFE = FEOperator(_res_y,_jac_y,X_diag,Y_diag)
  nls = GridapSolvers.NonlinearSolvers.NewtonSolver(ls_diag,verbose=i_am_main(ranks))
  yh0 = solve(nls,_opFE)
  psave(diag_dir*"/solT_$(t0)",yh0)


  #### PROGNOSTIC VARIABLES

  # equation for depth and velocity:
  mass(t,(dut,dpt),(v,r)) = ∫( (dut⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( (dpt*r)*meas_cf )dΩ

  res_p(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) = ∫( r*(F⋅grad_meas_cf + meas_cf*(∇⋅F) )  )dΩ

  res_u(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) = (  ∫( q*( (perp_matrix_cf⋅F) ⋅(metric_cf ⋅v))   )dΩ
                                + ∫( -τ*( (q-q0)/dt )*( (perp_matrix_cf⋅F) ⋅(metric_cf ⋅v))   )dΩ
                                + ∫( -τ*(u⋅∇(q))*( (perp_matrix_cf⋅F) ⋅(metric_cf ⋅v))   )dΩ
                                - ∫( Φ*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
                    )

  res_x(t,((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) = res_u(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) + res_p(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0))
  jac_x(t,((u,p),(q,F,Φ)),(du,dp),(v,r),(q0,F0,Φ0)) =  ∫( -τ*(du⋅∇(q))*( (perp_matrix_cf⋅F) ⋅(metric_cf ⋅v))   )dΩ
  jac_xt(t,((u,p),(q,F,Φ)),(dut,dpt),(v,r),(q0,F0,Φ0)) =  ∫( (dut⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( (dpt*r)*meas_cf )dΩ


  opT = TransientSemilinearFEOperator(mass,res_x,(jac_x,jac_xt),X_prog,Y_prog)
  opFE = FEOperator(res_y,jac_y,X_diag,Y_diag)
  opDAE = DAEFEOperator(opT,opFE,ls_diag)

  # transient parameters
  _dt = dx(nc(panel_model))*CFL/(p_fe*sqrt(gravity*_H_0))
  dt = floor(_dt, sigdigits=1)
  τ = dt/2

  # solve with SSP RK 3
  ode_solver = RungeKutta(ls_ode,ls_ode,dt,:EXRK_SSP_3_3)
  solT  = solve(ode_solver,opDAE,t0,_tF,xh0)

  ## iterate solution
  it = iterate(solT)

  unwrap_sw(it,ranks,solT,dir,_tF)


end

function unwrap_sw(it,ranks,solT,dir,tF,freq=10)
  sim_dir = dir*"/sim_data"
  final_dir = dir*"/final_solution"
  prog_dir = sim_dir*"/prognostics"
  diag_dir = sim_dir*"/diagnostics"

  counter = 1
  while !isnothing(it)
    data, state = it
    t, xh = data
    odeopcache = state[2][5][2]
    yh = odeopcache.diagnostics

    i_am_main(ranks) && println("t = ", t)

    if mod(counter,freq) == 0
      psave(prog_dir*"/solT_$t",xh)
      psave(diag_dir*"/solT_$t",yh)
    end

    if t >= tF - Gridap.ODEs.ε
      i_am_main(ranks) && println("Saving final solution")
      psave(final_dir*"/solT_$t",xh)
      # psave(final_dir*"/solT_diagnostics_$t",yh)
    end

    counter = counter + 1
    it = iterate(solT, state)
  end

end


function post_process(panel_model,p_fe::Int,dir::String,f::Function,return_vtk=false)

  # get the ranks to help with storing/saving solution
  ranks = get_ranks(panel_model)

  sim_dir = dir*"/sim_data"
  prog_dir = sim_dir*"/prognostics"
  diag_dir = sim_dir*"/diagnostics"

  vtk_dir = dir*"/vtk_data"
  (i_am_main(ranks) && !isdir(vtk_dir) ) && mkdir(vtk_dir)

  dir_casimirs = _dir*"/casimirs"
  (i_am_main(ranks) && !isdir(dir_casimirs)) && mkdir(dir_casimirs)

  lvl = nref(nc(panel_model))

  ## finite element solver
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  R = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1)
  H = TransientTrialFESpace(R)

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TransientTrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TransientTrialFESpace(V)

  X_prog = TransientMultiFieldFESpace([U,P]) # u, p
  Y_prog = MultiFieldFESpace([V,Q]) # u, p

  X_diag = TransientMultiFieldFESpace([H,U,P]) # q, F, Φ
  Y_diag = MultiFieldFESpace([R,V,Q]) # q, F, Φ

  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  cor_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  gravity = _g

  cell_geo_map = geo_map_func(Ω_panel)

  labels = ["uh","ph","qh","Fh","Phih","vort"]
  function make_vtk(t::Float64,xh,yh,cell_geo_map)
    uh,ph = xh
    qh,Fh,Φh = yh
    vort = qh*ph - cor_cf
    panel_cfs = [covarient_basis_cf⋅uh, ph, qh, Fh, Φh, vort]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,vtk_dir*"/solT_$t" * ".vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  function casimirs(xh,yh,dΩ)
    uh,ph = xh
    qh,Fh,Φh = yh
    vort = qh*ph - cor_cf

    ens = sum(∫( (qh*qh*xh[2])*meas_cf  )dΩ)
    energy = sum(∫( (0.5*xh[2]*( xh[1] ⋅(metric_cf⋅xh[1])) + 0.5*gravity*xh[2]*xh[2] )*meas_cf )dΩ)
    _mass = sum( ∫( xh[2]*meas_cf )dΩ  )
    _vort = sum( ∫( vort*meas_cf )dΩ  )

    return _mass, energy, ens, _vort
  end


  folders = readdir(prog_dir)
  dfolders = readdir(diag_dir)
  simName = "solT"

  ## casimirs to store
  ts = Vector{Float64}(undef,length(folders))
  Masss = Vector{Float64}(undef,length(folders))
  Energys = Vector{Float64}(undef,length(folders))
  Enstropys = Vector{Float64}(undef,length(folders))
  Vorts = Vector{Float64}(undef,length(folders))

  for (i,(f,g)) in enumerate(zip(folders,dfolders))
    t = parse(Float64,f[length(simName)+2:length(f)])

    x =  pload(joinpath(prog_dir,f),ranks)
    xh = FEFunction(X_prog,x)

    y =  pload(joinpath(diag_dir,f),ranks)
    yh = FEFunction(X_diag,y)

    i_am_main(ranks) && println("t = ", t)

    ts[i] = t
    Masss[i], Energys[i], Enstropys[i], Vorts[i] = casimirs(xh,yh,dΩ)

    return_vtk && make_vtk(t,xh,yh,cell_geo_map)

    if mod(i,10) == 0
      dxx = dx(panel_model)
      output = @strdict ts Masss Energys Enstropys Vorts dxx
      i_am_main(ranks) && safesave(datadir(dir_casimirs, ("casimirs.jld2")), output)
    end

  end

  dxx = dx(panel_model)
  output = @strdict ts Masss Energys Enstropys Vorts dxx
  i_am_main(ranks) && safesave(datadir(dir_casimirs, ("casimirs.jld2")), output)

  _make_pvd_distributed(vtk_dir,"solT",1)

end

function convergence_post_process(panel_model,p_fe::Int,dir::String)

  # get the ranks to help with storing/saving solution
  ranks = get_ranks(panel_model)

  initial_dir = dir*"/initial_solution"
  final_dir = dir*"/final_solution"

  lvl = nref(nc(panel_model))

  ## finite element solver
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ_error = Measure(Ω_panel,6*(p_fe+1))

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TransientTrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TransientTrialFESpace(V)

  X_prog = TransientMultiFieldFESpace([U,P]) # u, p
  Y_prog = MultiFieldFESpace([V,Q]) # u, p

  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  f_folders = readdir(final_dir)
  i_folders = readdir(initial_dir)
  simName = "solT"

  ## casimirs to store
  f = f_folders[1]
  g = i_folders[1]

  t = parse(Float64,f[length(simName)+2:length(f)])

  x =  pload(joinpath(final_dir,f),ranks)
  xh = FEFunction(X_prog,x)
  uh,ph = xh

  x0 =  pload(joinpath(initial_dir,g),ranks)
  xh0 = FEFunction(X_prog,x0)
  uh0,ph0 = xh0

  uh_proj = covarient_basis_cf ⋅ uh
  uh0_proj = covarient_basis_cf ⋅ uh0
  e_u = l2( (uh0_proj - uh_proj),meas_cf,dΩ_error)
  e_p = l2((ph0 - ph),meas_cf,dΩ_error)


  ### convergence output for DrWatson
  dir_convergence = dir*"/convergence"
  (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

  n = nc(panel_model)
  dxx = dx(panel_model)
  output = @strdict e_u e_p n dxx p_fe lvl t
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("shallow_water_nref$(lvl)_p$p_fe.jld2")), output)


end
################################################################################
#### Main run for transient solution
################################################################################
function main_transient(distribute,nprocs;
  restart=false,options="",n_ref_lvls=4,p_fe=1,CFL=0.1,ζ=0.0,return_vtk=true)

  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("transient_wave_equation")

  h = panel_to_cartesian(h₀(ζ))
  vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
  f = panel_to_cartesian(f₀(ζ))
  b = panel_to_cartesian(topography)

  ls_diag = CGSolver(JacobiLinearSolver();rtol=1-12,verbose=false,name="diagnostic_solver")
  ls_ode = CGSolver(JacobiLinearSolver();rtol=1-12,verbose=false,name="ode_solver")
  lss = (ls_ode,ls_diag)

  # models = get_models(ranks,nprocs,n_ref_lvls;threedims=false,octree=octree)
  # panel_model = models[1]
  omodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=n_ref_lvls)
  panel_model = omodel.parametric_dmodel

  _dir = datadir("TransientShallowWater_checkpointing")
  (i_am_main(ranks) && !isdir(_dir)) && mkdir(_dir)

  dir = _dir*"/sol_p$(p_fe)_nref$n_ref_lvls"
  (i_am_main(ranks) && !isdir(dir) && return_vtk) && mkdir(dir)

  transient_shallow_water_solver(panel_model,p_fe,dir,h,vX,f,b,CFL,lss,restart)
  convergence_post_process(panel_model,p_fe,dir)
  post_process(panel_model,p_fe,dir,f,return_vtk)

  i_am_main(ranks) && println("--DONE--")
  @test true
end


MPI.Init()
nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
# ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

with_mpi() do distribute
  main_transient(distribute,nprocs;restart=false,n_ref_lvls=2)
end
