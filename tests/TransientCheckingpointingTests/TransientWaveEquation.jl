"""
solve the linearised wave equation in
∂ₜu + ∇ᵧ(φ) = 0
∂ₜφ + ∇ᵧ⋅u = 0
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

## initial conditions
vecX(XYZ) = zero(XYZ)
depth(XYZ) = 1.0 + 0.01*exp(-5*((1-XYZ[1])^2+(0-XYZ[2])^2+(0-XYZ[3])^2))


function transient_wave_solver(panel_model::Union{<:DiscreteModel{2,2},<:GridapDistributed.DistributedDiscreteModel{2,2}},
  p_fe::Int,dir::String,h::Function,vX::Function,CFL=0.1,ls=LUSolver(),tF=2*π,restart=false)

  # get the ranks to help with storing/saving solution
  ranks = get_ranks(panel_model)

  sim_dir = dir*"/sim_data"
  (i_am_main(ranks) && !isdir(sim_dir) ) && mkdir(sim_dir)

  final_dir = dir*"/final_solution"
  (i_am_main(ranks) && !isdir(final_dir) ) && mkdir(final_dir)

  ## finite element solver
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TransientTrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TransientTrialFESpace(V)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  ## initial conditions
  function initial_condition()
    i_am_main(ranks) && println("initial condition")
    h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
    vec_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
    xh0 = interpolate([vec_contra_cf,h_cf],X)
    t = 0.0
    psave(sim_dir*"/solT_$(t)",xh0)
    return t,xh0
  end

  simName = "solT"
  t0,xh0 = (restart) ? load_last(ranks,X,sim_dir,simName) : initial_condition()

  ## transient weak form
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  mass(t, (dtu,dtp), (v,q)) = ∫( (v⋅ (metric_cf⋅ dtu) )*meas_cf )dΩ + ∫( (q*dtp)*meas_cf )dΩ
  res(t,(u,p),(v,q)) =  ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) )  )dΩ - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
  jac(t,(u,p),(du,dp),(v,q)) = res(t,(du,dp),(v,q))
  jac_t(t,(u,p),(dut,dpt),(v,q)) =  ∫( (dut⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( (dpt*q)*meas_cf )dΩ


  opT = TransientSemilinearFEOperator(mass, res,(jac,jac_t), X, Y, constant_mass=true)

  # transient parameters
  _dt = dx(panel_model)*CFL/p_fe
  dt = floor(_dt, sigdigits=1)

  # solve with SSP RK 3
  solver = RungeKutta(ls, ls, dt, :EXRK_SSP_3_3)
  solT = solve(solver, opT, t0, tF, xh0)

  ## iterate solution
  it = iterate(solT)

  unwrap(it,ranks,solT,dir,tF,10)


end


function post_process(panel_model,p_fe::Int,dir::String,return_vtk=false)

  # get the ranks to help with storing/saving solution
  ranks = get_ranks(panel_model)

  sim_dir = dir*"/sim_data"
  (i_am_main(ranks) && !isdir(sim_dir) ) && mkdir(sim_dir)

  vtk_dir = dir*"/vtk_data"
  (i_am_main(ranks) && !isdir(vtk_dir) ) && mkdir(vtk_dir)

  lvl = nref(nc(panel_model))

  ## finite element solver
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TransientTrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TransientTrialFESpace(V)

  Y = MultiFieldFESpace([V, Q])
  X = MultiFieldFESpace([U, P])

  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)


  cell_geo_map = geo_map_func(Ω_panel)

  labels = ["uh","ph"]
  function make_vtk(t::Float64,xh,cell_geo_map)
    uh,ph = xh
    panel_cfs = [covarient_basis_cf⋅uh, ph]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,vtk_dir*"/solT_$t" * ".vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  function casimirs(xh,dΩ)
    uh,ph = xh
    _m = sum( ∫(  meas_cf*ph )dΩ)
    _E = sum( ∫( 0.5*( uh⋅(metric_cf⋅ uh) + ph*ph)*meas_cf )dΩ  )
    _s_divu = sum(∫(   divergence(meas_cf*uh) )dΩ)
    _divu = sum(∫(  divergence(uh)  )dΩ)
    return _m, _E, _s_divu, _divu
  end


  folders = readdir(sim_dir)
  simName = "solT"

  ## casimirs to store
  ts = Vector{Float64}(undef,length(folders))
  ms = Vector{Float64}(undef,length(folders))
  Es = Vector{Float64}(undef,length(folders))
  s_divus = Vector{Float64}(undef,length(folders))
  divus = Vector{Float64}(undef,length(folders))

  for (i,f) in enumerate(folders)
    t = parse(Float64,f[length(simName)+2:length(f)])

    x =  pload(joinpath(sim_dir,f),ranks)
    xh = FEFunction(X,x)

    i_am_main(ranks) && println("t = ", t)

    ts[i] = t
    ms[i], Es[i], s_divus[i], divus[i] = casimirs(xh,dΩ)

    return_vtk && make_vtk(t,xh,cell_geo_map)

    if mod(i,10) == 0
      dxx = dx(panel_model)
      output = @strdict ts ms Es s_divus divus dxx
      i_am_main(ranks) && safesave(datadir(dir, ("wave_equation_nref$(lvl)_p$p_fe.jld2")), output)
    end

  end

  dxx = dx(panel_model)
  output = @strdict ts ms Es s_divus divus dxx
  i_am_main(ranks) && safesave(datadir(dir, ("wave_equation_nref$(lvl)_p$p_fe.jld2")), output)

  _make_pvd_distributed(vtk_dir,"solT",1)

end


################################################################################
#### Main run for transient solution
################################################################################
function main_transient(distribute,nprocs;
  restart=false,options="",n_ref_lvls=4,p_fe=1,CFL=0.1,tF=2*π,return_vtk=true)

  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("transient_wave_equation")

  h = panel_to_cartesian(depth)
  vX = panel_to_cartesian(tangent_vec(vecX))
  ls = LUSolver()

  _dir = datadir("TransientWaveEquation_checkpointing")
  (i_am_main(ranks) && !isdir(_dir)) && mkdir(_dir)

  # models = get_models(ranks,nprocs,n_ref_lvls;threedims=false,octree=octree)
  # panel_model = models[1]
  omodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=n_ref_lvls)
  panel_model =omodel.parametric_dmodel

  dir = _dir*"/sol_p$(p_fe)_nref$n_ref_lvls"
  (i_am_main(ranks) && !isdir(dir) && return_vtk) && mkdir(dir)


  GridapPETSc.Init(args=split(options))

  ls = CGSolver(JacobiLinearSolver();maxiter=2000,verbose=i_am_main(ranks))
  # ls = LUSolver()
  # ls = PETScLinearSolver()

  transient_wave_solver(panel_model,p_fe,dir,h,vX,CFL,ls,tF,restart)
  post_process(panel_model,p_fe,dir,return_vtk)

  GridapPETSc.Finalize()
  GridapPETSc.gridap_petsc_gc()

  i_am_main(ranks) && println("--DONE--")
  @test true
end


# MPI.Init()
# nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))

# with_mpi() do distribute
#   main_transient(distribute,nprocs;restart=false,n_ref_lvls=3)
# end
