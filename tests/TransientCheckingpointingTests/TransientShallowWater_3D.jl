"""
solve the non-linear shallow water equations in 3D
∂ₜu + q × F^† + ∇ᵧ(Φ) = 0
∂ₜφ + ∇ᵧ⋅F = 0
F = φu
Φ = 0.5(u⋅u) + gᵣφ
q = 1/φ( ∇ᵧ × u  + f )
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
include("../Geophysical/Williamson2Test_3D.jl")
include("../Geophysical/CurlConformingFESpacesFixes.jl")

inv_jacobian(p) = x -> inv(forward_jacobian_3D(p)(x))
contra_v_3D(vecX::Function,p::Int) = x -> inv_jacobian(p)(x) ⋅ vecX(p)(x)
contra_v_3D(vecX::Function) = p -> contra_v_3D(vecX,p)

transpose_jacobian(p) = x -> transpose(forward_jacobian_3D(p)(x))
inv_tranpose_jacobian(p) = x -> inv(transpose_jacobian(p)(x))
contravariant_basis_3D(p) = x -> inv_tranpose_jacobian(p)(x)

covar_v_3D(vecX::Function,p::Int) = x -> transpose_jacobian(p)(x) ⋅ vecX(p)(x)
covar_v_3D(vecX::Function) = p -> covar_v_3D(vecX,p)


function transient_shallow_water_solver_3D(
  panel_model::GridapDistributed.DistributedDiscreteModel{3,3},
  p_fe::Int,dir::String,
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
  dΩ = Measure(Ω_panel,4*(p_fe+1))

  tags = ["top_boundary", "bottom_boundary"]
  Γ = BoundaryTriangulation(panel_model,tags=tags)
  dΓ = Measure(Γ,4*(p_fe+1))
  nΓ = get_normal_vector(Γ)

  R = TestFESpace(Ω_panel, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl,dirichlet_tags=tags)
  H = TrialFESpace(R,VectorValue(0.0,0.0,0.0))

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv,dirichlet_tags=tags)
  U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))

  X_prog = MultiFieldFESpace([U,P]) # u, p
  Y_prog = MultiFieldFESpace([V,Q]) # u, p

  X_diag = MultiFieldFESpace([H,U,P]) # q, F, Φ
  Y_diag = MultiFieldFESpace([R,V,Q]) # q, F, Φ

  ## initial conditions
  function initial_condition()
    i_am_main(ranks) && println("initial condition")

    u_contra_cf = panelwise_cellfield(contra_v_3D(u_vec_3D),Ω_panel,panel_ids)
    u_contra_h = interpolate(u_contra_cf,U)

    h_cf = panelwise_cellfield(h_3D,Ω_panel,panel_ids)
    b_cf = panelwise_cellfield(topography,Ω_panel,panel_ids)
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
  grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)
  inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  contravariant_basis_cf = panelwise_cellfield(contravariant_basis_3D,Ω_panel,panel_ids)
  jac_cf = panelwise_cellfield(forward_jacobian,Ω_panel,panel_ids)
  area_meas_cf = Operation(norm)(jac_cf⋅(inv_metric_cf ⋅nΓ) )

  gravity = _g
  f_cov_cf = panelwise_cellfield(covar_v_3D(f_vec_3D),Ω_panel,panel_ids)
  n_cov = CellField(x->VectorValue(1,0,0),Ω_panel) ## same as nΓ
  b_cf = panelwise_cellfield(topography,Ω_panel,panel_ids) # topography

  #### DIAGNOSTIC VARIABLES
  # vorticity
  resq(((u,p),(q,F,Φ)),(w,v,ψ)) = ( ∫( p*(q⋅(inv_metric_cf⋅w))*meas_cf )dΩ
                                  - ∫( u⋅( metric_cf⋅ curl(w) )  )dΩ
                                  + ∫( (( w × (metric_cf⋅ u) )⋅nΓ)*area_meas_cf   )dΓ
                                  - ∫( (f_cov_cf⋅(inv_metric_cf ⋅ w))*meas_cf )dΩ
                                  )

  # mass flux
  resF(((u,p),(q,F,Φ)),(w,v,ψ)) = ∫( (F⋅ (metric_cf⋅v))*meas_cf )dΩ - ∫( p*(u⋅(metric_cf⋅v))*meas_cf   )dΩ

  # Bernoulli potential
  resΦ(((u,p),(q,F,Φ)),(w,v,ψ)) = ∫( Φ*ψ*meas_cf  )dΩ - ∫( gravity*(p+b_cf)*ψ*meas_cf  )dΩ - ∫( 0.5*( u ⋅(metric_cf⋅u) )ψ*meas_cf  )dΩ

  res_y(t,((u,p),(q,F,Φ)),(w,v,ψ)) = resq(((u,p),(q,F,Φ)),(w,v,ψ)) + resF(((u,p),(q,F,Φ)),(w,v,ψ)) + resΦ(((u,p),(q,F,Φ)),(w,v,ψ))
  jac_y(t,((u,p),(q,F,Φ)),(dq,dF,dΦ),(w,v,ψ)) = ∫( p*(dq⋅(inv_metric_cf⋅w))*meas_cf )dΩ  + ∫( (dF⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( dΦ*ψ*meas_cf  )dΩ

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

  res_u(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) = (
                                  ∫( ( q × (metric_cf⋅F) )⋅(metric_cf⋅v) )dΩ
                                + ∫( -τ*( ( (q-q0)/dt ) × (metric_cf⋅F) )⋅(metric_cf⋅v) )dΩ
                                + ∫( τ*(  ( u× curl(q)) × (metric_cf⋅F)  )⋅(metric_cf⋅v)   )dΩ
                                - ∫( Φ*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
                    )

  res_x(t,((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) = res_u(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) + res_p(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0))
  jac_x(t,((u,p),(q,F,Φ)),(du,dp),(v,r),(q0,F0,Φ0)) = + ∫( τ*(  ( du× curl(q)) × (metric_cf⋅F)  )⋅(metric_cf⋅v)   )dΩ
  jac_xt(t,((u,p),(q,F,Φ)),(dut,dpt),(v,r),(q0,F0,Φ0)) =  ∫( (dut⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( (dpt*r)*meas_cf )dΩ


  opT = TransientSemilinearFEOperator(mass,res_x,(jac_x,jac_xt),X_prog,Y_prog)
  opFE = FEOperator(res_y,jac_y,X_diag,Y_diag)
  opDAE = DAEFEOperator(opT,opFE,ls_diag)

  # transient parameters
  dxx_horizontal = dx_horizontal(panel_model)
  _dt = dxx_horizontal*CFL/(p_fe*sqrt(gravity*_H_0))
  dt = floor(_dt, sigdigits=1)
  i_am_main(ranks) && println("dt = $dt")
  τ = dt/2

  # solve with SSP RK 3
  ode_solver = RungeKutta(ls_ode,ls_ode,dt,:EXRK_SSP_3_3)
  solT  = solve(ode_solver,opDAE,t0,_tF,xh0)

  ## iterate solution
  it = iterate(solT)

  unwrap_sw(it,ranks,solT,dir,_tF)


end

function unwrap_sw(it,ranks,solT,dir,tF,freq=25)
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


function post_process(panel_model,p_fe::Int,dir::String,return_vtk=true)

  # get the ranks to help with storing/saving solution
  ranks = get_ranks(panel_model)

  sim_dir = dir*"/sim_data"
  prog_dir = sim_dir*"/prognostics"
  diag_dir = sim_dir*"/diagnostics"

  vtk_dir = dir*"/vtk_data"
  (i_am_main(ranks) && !isdir(vtk_dir) ) && mkdir(vtk_dir)

  vtk_latlon_dir = dir*"/latlon_data"
  (i_am_main(ranks) && !isdir(vtk_latlon_dir) ) && mkdir(vtk_latlon_dir)

  dir_casimirs = dir*"/casimirs"
  (i_am_main(ranks) && !isdir(dir_casimirs)) && mkdir(dir_casimirs)

  # ensure no MPI task tries to generate the file before the main MPI task has
  # created the folder
  PartitionedArrays.barrier(ranks)

  ## finite element solver
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,4*(p_fe+1))

  tags = ["top_boundary", "bottom_boundary"]
  R = TestFESpace(panel_model, ReferenceFE(nedelec,Float64,p_fe);conformity=:Hcurl,dirichlet_tags=tags)
  H = TrialFESpace(R,VectorValue(0.0,0.0,0.0))

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TransientTrialFESpace(Q)

  V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv,dirichlet_tags=tags)
  U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))

  X_prog = TransientMultiFieldFESpace([U,P]) # u, p
  Y_prog = MultiFieldFESpace([V,Q]) # u, p

  X_diag = TransientMultiFieldFESpace([H,U,P]) # q, F, Φ
  Y_diag = MultiFieldFESpace([R,V,Q]) # q, F, Φ

  inv_metric_cf = panelwise_cellfield(inv_metric,Ω_panel,panel_ids)
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  f_cov_cf = panelwise_cellfield(covar_v_3D(f_vec_3D),Ω_panel,panel_ids)
  gravity = _g

  cell_geo_map = geo_map_func(Ω_panel)
  latlon_geo_map = latlon_geo_map_func(Ω_panel)
  labels = ["uh","ph","qh","Fh","Phih","vort"]
  function make_vtk(t::Float64,xh,yh,cell_geo_map,latlon_geo_map)
    uh,ph = xh
    qh,Fh,Φh = yh
    # vort = qh*ph - f_cov_cf
    panel_cfs = [covarient_basis_cf⋅uh, ph,
                 covarient_basis_cf ⋅ (inv_metric_cf ⋅ qh ),
                 Fh, Φh,
                #  covarient_basis_cf ⋅ (inv_metric_cf ⋅ vort )
                 ]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,vtk_dir*"/solT_$t" * ".vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)
    writevtk(Ω_panel,vtk_latlon_dir*"/solT_$t" * ".vtu", cellfields=cellfields,append=false,geo_map=latlon_geo_map)
  end

  function casimirs(xh,yh,dΩ)
    uh,ph = xh
    qh,Fh,Φh = yh
    vort = 0.0#qh*ph - f_cov_cf

    ens = 0.0#sum(∫( ( (qh⋅qh)*xh[2])*meas_cf  )dΩ)
    energy = sum(∫( (0.5*xh[2]*( xh[1] ⋅(metric_cf⋅xh[1])) + 0.5*gravity*xh[2]*xh[2] )*meas_cf )dΩ)
    _mass = sum( ∫( xh[2]*meas_cf )dΩ  )
    _vort = 0.0#sum( ∫( vort*meas_cf )dΩ  )

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

    return_vtk && make_vtk(t,xh,yh,cell_geo_map,latlon_geo_map)

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
  _make_pvd_distributed(vtk_latlon_dir,"solT",1)
end

################################################################################
#### Main run for transient solution
################################################################################
function main_transient(distribute,nprocs;
  restart=false,options="",n_ref_lvls=4,p_fe=1,CFL=0.1,return_vtk=true)

  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("transient_wave_equation_3D")

  options = """
  -cg_ksp_type cg
  -cg_ksp_converged_reason
  -cg_ksp_monitor
  -cg_ksp_rtol 1.0e-10
  -cg_pc_type gamg
  """

  o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
        num_horizontal_uniform_refinements=n_ref_lvls,
        num_vertical_uniform_refinements=0)
  panel_model = o3model.parametric_dmodel

  # ls_diag = CGSolver(JacobiLinearSolver();rtol=1-8,atol=1e-10,verbose=i_am_main(ranks),name="diagnostic_solver")
  # ls_ode = CGSolver(JacobiLinearSolver();rtol=1-8,atol=1e-10,verbose=i_am_main(ranks),name="ode_solver")

  GridapPETSc.Init(args=split(options))
  ls_ode = PETScLinearSolver(petsc_cg_amg_setup)
  ls_diag = PETScLinearSolver(petsc_cg_amg_setup)

  lss = (ls_ode,ls_diag)

  _dir = datadir("TransientShallowWater_3D")
  (i_am_main(ranks) && !isdir(_dir)) && mkdir(_dir)

  dir = _dir*"/sol_p$(p_fe)_nref$n_ref_lvls"
  (i_am_main(ranks) && !isdir(dir) && return_vtk) && mkdir(dir)

  ## if restarted, post process the existing files
  restart && post_process(panel_model,p_fe,dir,return_vtk)
  i_am_main(ranks) && println("finished post processing existing data")

  transient_shallow_water_solver_3D(panel_model,p_fe,dir,CFL,lss,restart)
  post_process(panel_model,p_fe,dir,return_vtk)

  GridapPETSc.Finalize()
  GridapPETSc.gridap_petsc_gc()

  i_am_main(ranks) && println("--DONE--")
  @test true
end

function main_visualise(distribute,nprocs;
  n_ref_lvls=4,p_fe=1,return_vtk=true)

  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("Transient SW 3D visualise")

  o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
  num_horizontal_uniform_refinements=n_ref_lvls,
  num_vertical_uniform_refinements=0)
  panel_model = o3model.parametric_dmodel

  _dir = datadir("TransientShallowWater_3D")

  dir = _dir*"/sol_p$(p_fe)_nref$n_ref_lvls"

  post_process(panel_model,p_fe,dir,return_vtk)

  i_am_main(ranks) && println("--DONE--")
  @test true
end


function petsc_cg_amg_setup(ksp)
  rtol = GridapPETSc.PETSC.PETSC_DEFAULT
  atol = GridapPETSc.PETSC.PETSC_DEFAULT
  dtol = GridapPETSc.PETSC.PETSC_DEFAULT
  maxits = PetscInt(1000)

  @check_error_code GridapPETSc.PETSC.KSPSetOptionsPrefix(ksp[],"cg_")
  @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])
  # @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[],GridapPETSc.PETSC.KSPGMRES)

  pc = Ref{GridapPETSc.PETSC.PC}()
  @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[],pc)
  @check_error_code GridapPETSc.PETSC.PCSetType(pc[],GridapPETSc.PETSC.PCGAMG)
  @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[], rtol, atol, dtol, maxits)
  @check_error_code GridapPETSc.PETSC.KSPView(ksp[],C_NULL)
end


# MPI.Init()
# nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
# # ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))
# # n_ref_lvls = 3
# with_mpi() do distribute
#   main_transient(distribute,nprocs;restart=false,n_ref_lvls=3)
# end
