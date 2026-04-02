 """
solve the transient linear Boussineq equations in 3D
∂ₜu + f(̂n×u) + ∇ᵧ(φ) - bn̂ = 0.0
∂ₜφ + c² ∇ᵧ⋅u = 0.0
∂ₜb + N² u⋅̂n = 0.0
"""

using MPI
using PartitionedArrays

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Test
using GridapPETSc

import GridapGeosciences.Helpers: RADIUS, THICKNESS


include("helpers.jl")
include("../convergence_tools.jl")

a_e = 6.37e6/125 # m
d = 5000 #m
Lz = 20e3 #m
R = a_e # m radius
u_0 = 20 #m/s
Ωr = 7.292e-5 #1/s
c = 343 #m/s speed of sound
N = 0.01 #1/s bouyancy frequency
ztop = 10e3 #m
dΘ = 1 #K
TF = (3600*24)*10 # s

LH = a_e # m
LV = ztop/THICKNESS
c2 = c/LH # 1/s
τ = 1/c2#1/Ωr # s

_R = R/LH
_ztop = ztop/LV
_d = d/LV
_Lz = Lz/LV
_u_0 = u_0*τ/LH
_Ωr = Ωr*τ
_c = c*τ/LH
_N = N*τ
tF = TF/τ

p0(xyz) = 0.0

function b0(xyz)
  x,y,z = xyz
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr

  θc = 2*π/3
  ϕc = 0.0

  k = sqrt(x^2 + y^2 + z^2) - _R

  r = _R*acos( sin(ϕc)*sin(ϕ) + cos(ϕc)*cos(ϕ)*cos(θ-θc)    )
  s = _d^2/(_d^2 + r^2)
  b = dΘ*s*sin( 2*π*k/_Lz  )
  b
end

function u0(xyz)
  x,y,z = xyz
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr

  # u = _u_0*cos(ϕ) #
  # v = 0.0#
  u = -_u_0*y/_R
  v = _u_0*x/_R

  VectorValue(u,v,0.0)
end

function omega(xyz)
  x,y,z = xyz
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  2*_Ωr*sin(ϕ)
end


function transient_linear_boussinesq_solver(
  panel_model::GridapDistributed.DistributedDiscreteModel{3,3},
  p_fe::Int,dir::String,h::Function,vX::Function,f::Function,b::Function,
  ls=LUSolver(),CFL=0.1,restart=false)

  ranks = get_ranks(panel_model)


  lvl_h = nref(nc_horizontal(panel_model))
  lvl_v = nref(nc_vertical(panel_model))
  i_am_main(ranks) && println("nref_h = $lvl_h; nref_v = $lvl_v; p_fe = $p_fe")

  sim_dir = dir*"/sim_data"
  (i_am_main(ranks) && !isdir(sim_dir) ) && mkdir(sim_dir)

  final_dir = dir*"/final_solution"
  (i_am_main(ranks) && !isdir(final_dir) ) && mkdir(final_dir)

  initial_dir = dir*"/initial_solution"
  (i_am_main(ranks) && !isdir(initial_dir) ) && mkdir(initial_dir)

  # ensure no MPI task tries to generate the file before the main MPI task has
  # created the folder
  PartitionedArrays.barrier(ranks)

  ## finite element solver
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,4*(p_fe+1))

  tags = ["bottom_boundary",  "top_boundary"]

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv,dirichlet_tags=tags)
  U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))

  W = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  B = TrialFESpace(W)

  Y = MultiFieldFESpace([V, Q, W])
  X = MultiFieldFESpace([U, P, B])

  ## initial conditions
  function initial_condition()
    i_am_main(ranks) && println("initial condition")

    h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
    u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
    b_cf = panelwise_cellfield(b,Ω_panel,panel_ids)

    xh0 = interpolate([u_contra_cf,h_cf,b_cf],X)
    t = 0.0
    psave(sim_dir*"/solT_$(t)",xh0)
    psave(initial_dir*"/solT_$(t)",xh0)
    return t,xh0
  end

  simName = "solT"
  t0,xh0 = (restart) ? load_last(ranks,X,sim_dir,simName) : initial_condition()

  omega_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  g_star_cf = panelwise_cellfield(g_star,Ω_panel,panel_ids)

  # weak forms
  detg_cf = panelwise_cellfield(detg,Ω_panel,panel_ids)
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)

  #### Velocity
  Aperp = [0 0 0
            0 0 -1
            0 1 0]
  Rperp = TensorValue(Aperp)
  Rperp_cf = CellField(Rperp,Ω_panel)

  mass(t, (dtu,dtp,dtb), (v,q,r)) = ( ∫( (v⋅ (metric_cf⋅ dtu) )*meas_cf )dΩ
                                    + ∫( (q*dtp)*meas_cf )dΩ
                                    + ∫( (r*dtb)*meas_cf )dΩ )
  #### Velocity
  resu(t,(u,p,b),(v,q,r)) = ( ∫( ( omega_cf*( (Rperp_cf⋅ u)⋅v))*detg_cf )dΩ
                              - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
                              - ∫( b*(g_star_cf⋅v )*meas_cf )dΩ )

  #### Pressure
  resp(t,(u,p,b),(v,q,r)) = ∫( _c^2*( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) ) ) )dΩ

  #### Bouyancy
  n_cf = CellField(VectorValue(1,0,0),Ω_panel)
  resb(t,(u,p,b),(v,q,r)) =  ∫( _N^2*( r*(g_star_cf⋅u)*meas_cf)   )dΩ
                                                  # ∫( _N^2*( r*( n_cf⋅(metric_cf⋅u))*meas_cf)   )dΩ

  res(t,(u,p,b),(v,q,r)) = resu(t,(u,p,b),(v,q,r)) + resp(t,(u,p,b),(v,q,r)) + resb(t,(u,p,b),(v,q,r))
  jac(t,(u,p,b),(du,dp,db),(v,q,r)) = resu(t,(du,dp,db),(v,q,r)) + resp(t,(du,dp,db),(v,q,r)) + resb(t,(du,dp,db),(v,q,r))
  jac_t(t,(u,p,b),(dut,dpt,dbt),(v,q,r)) =  mass(t, (dut,dpt,dbt), (v,q,r))

  opT = TransientSemilinearFEOperator(mass, res, (jac,jac_t), X, Y; constant_mass=true)

  # transient parameters
  dxx_horizontal = dx_horizontal(panel_model)
  _dt = dxx_horizontal*CFL/_c
  dt = floor(_dt, sigdigits=1)

  dxx_vertical = dx_vertical(panel_model)
  dxx_horizontal/dxx_vertical

  # solve with SSP RK 3
  nls = GridapSolvers.NonlinearSolvers.NewtonSolver(ls;verbose=i_am_main(ranks))
  # solver =  BackwardEuler(nls, dt)
  # solver = RungeKutta(ls,ls, dt,:EXRK_SSP_3_3)
  solver = ThetaMethod(nls, dt, 0.5)
  solT = solve(solver, opT, t0, tF, xh0)

  ## iterate solution
  it = iterate(solT)
  unwrap_lbous(it,ranks,solT,dir,tF)
end


function unwrap_lbous(it,ranks,solT,dir,tF,freq=25)
  sim_dir = dir*"/sim_data"
  final_dir = dir*"/final_solution"

  counter = 1
  while !isnothing(it)
    data, state = it
    t, xh = data

    i_am_main(ranks) && println("t = ", t)

    if mod(counter,freq) == 0
      psave(sim_dir*"/solT_$t",xh)
    end

    if t >= tF - Gridap.ODEs.ε
      i_am_main(ranks) && println("Saving final solution")
      psave(final_dir*"/solT_$t",xh)
    end

    counter = counter + 1
    it = iterate(solT, state)
  end

end

function post_process(panel_model,p_fe::Int,dir::String,return_vtk=true)

  # get the ranks to help with storing/saving solution
  ranks = get_ranks(panel_model)

  sim_dir = dir*"/sim_data"

  vtk_dir = dir*"/vtk_data"
  (i_am_main(ranks) && !isdir(vtk_dir) ) && mkdir(vtk_dir)

  vtk_latlon_dir = dir*"/vtk_latlon_data"
  (i_am_main(ranks) && !isdir(vtk_latlon_dir) ) && mkdir(vtk_latlon_dir)

  # ensure no MPI task tries to generate the file before the main MPI task has
  # created the folder
  PartitionedArrays.barrier(ranks)

  ## finite element solver
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)

  tags = ["bottom_boundary",  "top_boundary"]

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv,dirichlet_tags=tags)
  U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))

  W = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  B = TrialFESpace(W)

  Y = MultiFieldFESpace([V, Q, W])
  X = MultiFieldFESpace([U, P, B])

  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  cell_geo_map = geo_map_func(Ω_panel)
  latlon_cell_geo_map = latlon_geo_map_func(Ω_panel)
  # owned_panel_ids = get_owned_panel_ids(panel_model)

  # latlon_cell_geo_map = map(owned_panel_ids) do pid
  #   cell_geo_map = lazy_map(p -> ForwardMap(p), pid)
  #   fi = lazy_map(p->Cartesian2SphereicalMap3D(),pid)
  #   return lazy_map(∘, fi, cell_geo_map)
  # end

  labels = ["uh","ph", "bh"]
  function make_vtk(t::Float64,xh,cell_geo_map,latlon_cell_geo_map)
    uh,ph,bh = xh
    panel_cfs = [covarient_basis_cf⋅uh, ph, bh]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,vtk_dir*"/solT_$t",cellfields=cellfields,append=false,geo_map=cell_geo_map)
    writevtk(Ω_panel,vtk_latlon_dir*"/solT_$t",cellfields=cellfields,append=false,geo_map=latlon_cell_geo_map)
  end

  folders = readdir(sim_dir)
  simName = "solT"

  # ts = Vector{Float64}(undef,length(folders))
  for (i,f) in enumerate(folders)
    t = parse(Float64,f[length(simName)+2:length(f)])

    x =  pload(joinpath(sim_dir,f),ranks)
    xh = FEFunction(X,x)

    i_am_main(ranks) && println("t = ", t)

    # ts[i] = t

    return_vtk && make_vtk(t,xh,cell_geo_map,latlon_cell_geo_map)
  end

  # _make_pvd_distributed(vtk_dir,"solT",1)
  # _make_pvd_distributed(vtk_latlon_dir,"solT",1)

end

################################################################################
#### Main run for transient solution
################################################################################
function main_transient(distribute,nprocs;
  restart=false,n_ref_lvls=4,p_fe=1,CFL=0.1,return_vtk=true)
  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("transient_linear_boussinesq")

  h = panel_to_cartesian(p0)
  vX = panel_to_cartesian(tangent_vec(u0))
  f = panel_to_cartesian(omega)
  b = panel_to_cartesian(b0)

  options = """
      -g_ksp_type gmres
      -g_ksp_converged_reason
      -g_ksp_monitor
      """

  GridapPETSc.Init(args=split(options))
  ls = PETScLinearSolver(petsc_gmres_amg_setup)

  # ls = GMRESSolver(10;Pr=JacobiLinearSolver(),maxiter=1000,verbose=i_am_main(ranks))

  o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
          num_horizontal_uniform_refinements=n_ref_lvls,
          num_vertical_uniform_refinements=n_ref_lvls)

  panel_model = o3model.parametric_dmodel

  _dir = datadir("TransientLinearisedBoussinesq")
  (i_am_main(ranks) && !isdir(_dir)) && mkdir(_dir)

  dir = _dir*"/sol_p$(p_fe)_nref_h$(n_ref_lvls)_nref_v$(n_ref_lvls)"
  (i_am_main(ranks) && !isdir(dir) && return_vtk) && mkdir(dir)

  ## if restarted, post process the existing files
  restart && post_process(panel_model,p_fe,dir,return_vtk)
  i_am_main(ranks) && println("finished post processing existing data")

  transient_linear_boussinesq_solver(panel_model,p_fe,dir,h,vX,f,b,ls,CFL,restart)
  post_process(panel_model,p_fe,dir,return_vtk)

  GridapPETSc.Finalize()
  GridapPETSc.gridap_petsc_gc()

  i_am_main(ranks) && println("--DONE--")

end

function petsc_gmres_amg_setup(ksp)
  rtol = GridapPETSc.PETSC.PETSC_DEFAULT
  atol = GridapPETSc.PETSC.PETSC_DEFAULT
  dtol = GridapPETSc.PETSC.PETSC_DEFAULT
  maxits = PetscInt(1000)

  @check_error_code GridapPETSc.PETSC.KSPSetOptionsPrefix(ksp[],"g_")
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
# ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

# with_mpi() do distribute
#   main_transient(distribute,nprocs;restart=true,n_ref_lvls=3)
# end
