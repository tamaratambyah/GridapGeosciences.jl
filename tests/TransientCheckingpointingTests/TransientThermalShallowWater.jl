"""
solve the thermal shallow water equations
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
include("../Geophysical/ThermogeostrophicBalanceTest.jl")

function my_mean( Bu_n::SkeletonPair)
  plus  = ( Bu_n.plus)
  minus = ( Bu_n.minus)
  0.5*( plus - minus  )
end

# function upwinding_sign(Fn)
#   c = 0.0

#   if Fn < 0.0
#     c = -0.5
#   elseif Fn > 0.0
#     c = 0.5
#   end
#   return c
# end


function transient_tsw_solver(panel_model::Union{<:DiscreteModel{2,2},<:GridapDistributed.DistributedDiscreteModel{2,2}},
  p_fe::Int,dir::String,h::Function,vX::Function,f::Function,B::Function,
  ε=1e-4,soft=false,
  CFL=0.1,lss=(LUSolver(),LUSolver()),restart=false
  )

  # upwinding function
  function upwinding_sign(Fn)
    c = 0.0

    if Fn < -ε
      c = -0.5
    elseif Fn > ε
      c = 0.5
    end

    if soft
      c = 0.5*Fn/(sqrt(Fn^2 + (ε)^2 ) )
    end

    return c

  end

  ls_ode, ls_diag = lss

  das = FullyAssembledRows()

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

  # ensure no MPI task tries to generate the file before the main MPI task has
  # created the folder
  PartitionedArrays.barrier(ranks)

  ## finite element solver
  degree = 2*(p_fe+1)
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(das,panel_model)
  dΩ = Measure(Ω_panel,degree)

  Λ = SkeletonTriangulation(das,panel_model)
  dΛ = Measure(Λ,degree)
  n_Λ = get_normal_vector(Λ)

  R = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1)
  H = TransientTrialFESpace(R)

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TransientTrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TransientTrialFESpace(V)

  X_prog = TransientMultiFieldFESpace([U,P,P]) # u, p, B
  Y_prog = MultiFieldFESpace([V,Q,Q]) # u, p, B

  X_diag = TransientMultiFieldFESpace([H,U,P,P]) # q, F, Φ, b
  Y_diag = MultiFieldFESpace([R,V,Q,Q]) # q, F, Φ, b

  ## initial conditions
  function initial_condition()
    i_am_main(ranks) && println("initial condition")

    covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
    u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
    u_contra_h = interpolate(u_contra_cf,U)
    u_proj_h = covarient_basis_cf ⋅ u_contra_h

    h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
    B_cf = panelwise_cellfield(B,Ω_panel,panel_ids)
    h_h = interpolate(h_cf,P)
    B_h = interpolate(B_cf,P)

    xh0 = interpolate_everywhere([u_contra_h,h_h,B_h],X_prog(0.0))
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
  meas_cf_skel = panelwise_cellfield(sqrtg,Λ)

  #### DIAGNOSTIC VARIABLES
  #### T = 0.5p
  assem_diag = SparseMatrixAssembler(X_diag,Y_diag,das)

  # vorticity
  perp_matrix_cf = panelwise_cellfield(perp_matrix,Ω_panel,panel_ids)
  resq(((u,p,B),(q,F,Φ,b)),(w,v,ψ,r)) = ∫( q*p*w*meas_cf  )dΩ - ∫( cor_cf*w*meas_cf  )dΩ - ∫( (perp_matrix_cf⋅u)⋅∇(w)  )dΩ

  # mass flux
  resF(((u,p,B),(q,F,Φ,b)),(w,v,ψ,r)) = ∫( (F⋅ (metric_cf⋅v))*meas_cf )dΩ - ∫( p*(u⋅(metric_cf⋅v))*meas_cf   )dΩ

  # Bernoulli potential
  resΦ(((u,p,B),(q,F,Φ,b)),(w,v,ψ,r)) = ∫( Φ*ψ*meas_cf  )dΩ - ∫( 0.5*B*ψ*meas_cf  )dΩ - ∫( 0.5*( u ⋅(metric_cf⋅u) )ψ*meas_cf  )dΩ

  # Bouyancy
  resb(((u,p,B),(q,F,Φ,b)),(w,v,ψ,r)) = ∫( b*p*r*meas_cf  )dΩ - ∫( B*r*meas_cf  )dΩ


  res_y(t,((u,p,B),(q,F,Φ,b)),(w,v,ψ,r)) = (
      resq(((u,p,B),(q,F,Φ,b)),(w,v,ψ,r))
    + resF(((u,p,B),(q,F,Φ,b)),(w,v,ψ,r))
    + resΦ(((u,p,B),(q,F,Φ,b)),(w,v,ψ,r))
    + resT(((u,p,B),(q,F,Φ,b)),(w,v,ψ,r))
    + resb(((u,p,B),(q,F,Φ,b)),(w,v,ψ,r))
  )
  jac_y(t,((u,p,B),(q,F,Φ,b)),(dq,dF,dΦ,db),(w,v,ψ,r)) = (
      ∫( dq*p*w*meas_cf  )dΩ
    + ∫( (dF⋅ (metric_cf⋅v))*meas_cf )dΩ
    + ∫( dΦ*ψ*meas_cf  )dΩ
    + ∫( db*p*r*meas_cf  )dΩ
  )

  _res_y((q,F,Φ,b),(w,v,ψ,r))  = res_y(0.0,(xh0,(q,F,Φ,b)),(w,v,ψ,r))
  _jac_y((q,F,Φ,b),(dq,dF,dΦ,db),(w,v,ψ,r)) = jac_y(0.0,(xh0,(q,F,Φ,b)),(dq,dF,dΦ,db),(w,v,ψ,r))
  _opFE = FEOperator(_res_y,_jac_y,X_diag,Y_diag,assem_diag)
  nls = GridapSolvers.NonlinearSolvers.NewtonSolver(ls_diag,verbose=i_am_main(ranks))
  yh0 = solve(nls,_opFE)
  psave(diag_dir*"/solT_$(t0)",yh0)


  #### PROGNOSTIC VARIABLES
  assem_prog = SparseMatrixAssembler(X_prog,Y_prog,das)

  # equation for depth and velocity:
  mass(t,(dut,dpt,dBt),(v,r,w)) = (
      ∫( (dut⋅ (metric_cf⋅v))*meas_cf )dΩ
    + ∫( (dpt*r)*meas_cf )dΩ
    + ∫( (dBt*w)*meas_cf )dΩ
  )

  res_p(((u,p,B),(q,F,Φ,b)),(v,r,w),(q0,F0,Φ0,b0)) = ∫( r*(F⋅grad_meas_cf + meas_cf*(∇⋅F) )  )dΩ

  res_u(((u,p,B),(q,F,Φ,b)),(v,r,w),(q0,F0,Φ0,b0)) = (
            ∫( q*( (perp_matrix_cf⋅F) ⋅(metric_cf ⋅v))   )dΩ
          - ∫( Φ*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
          + ∫( 0.5*(b*(∇(0.5*p)⋅v) )*meas_cf )dΩ
          + ∫( -0.5*(b*(0.5*p))*(v⋅grad_meas_cf + meas_cf*(∇⋅v) )  )dΩ
          + ∫( -0.5*((0.5*p)*(∇(b)⋅v) )*meas_cf )dΩ
      )

  u_s1(((u,p,B),(q,F,Φ,b)),(v,r,w)) = (
      ∫( -0.5*my_mean((v*b)⋅n_Λ)*jump(0.5*p)*meas_cf_skel.plus   )dΛ
    + ∫( 0.5*my_mean((v*(0.5*p))⋅n_Λ)*jump(b)*meas_cf_skel.plus   )dΛ
  )

  u_s2(((u,p,B),(q,F,Φ,b)),(v,r,w)) = ∫( -0.5*( (upwinding_sign∘((F⋅ n_Λ).plus))*(v⋅n_Λ).plus )*jump(b)*jump(0.5*p)*meas_cf_skel.plus   )dΛ


  res_B(((u,p,B),(q,F,Φ,b)),(v,r,w),(q0,F0,Φ0,b0)) = (
      ∫( -0.5*(b*(∇(w)⋅F) )*meas_cf )dΩ
    + ∫( 0.5*(b*w)*(F⋅grad_meas_cf + meas_cf*(∇⋅F) )  )dΩ
    + ∫( 0.5*(w*(∇(b)⋅F) )*meas_cf )dΩ
  )



  B_s1(((u,p,B),(q,F,Φ,b)),(v,r,w)) = (
      ∫( 0.5*my_mean((F*b)⋅n_Λ)*jump(w)*meas_cf_skel.plus   )dΛ
    + ∫( 0.5*my_mean((F*w)⋅n_Λ)*jump(b)*meas_cf_skel.plus   )dΛ
  )

  B_s2(((u,p,B),(q,F,Φ,b)),(v,r,w)) = ∫( 0.5*( (upwinding_sign∘((F⋅ n_Λ).plus))*(F⋅n_Λ).plus )*jump(b)*jump(w)*meas_cf_skel.plus   )dΛ


  res_x(t,((u,p,B),(q,F,Φ,b)),(v,r,w),(q0,F0,Φ0,b0)) = (
      res_u(((u,p,B),(q,F,Φ,b)),(v,r,w),(q0,F0,Φ0,b0))
    + u_s1(((u,p,B),(q,F,Φ,b)),(v,r,w))
    + u_s2(((u,p,B),(q,F,Φ,b)),(v,r,w))
    + res_p(((u,p,B),(q,F,Φ,b)),(v,r,w),(q0,F0,Φ0,b0))
    + res_B(((u,p,B),(q,F,Φ,b)),(v,r,w),(q0,F0,Φ0,b0))
    + B_s1(((u,p,B),(q,F,Φ,b)),(v,r,w))
    + B_s2(((u,p,B),(q,F,Φ,b)),(v,r,w))
  )
  jac_xt(t,((u,p,B),(q,F,Φ,b)),(dut,dpt,dBt),(v,r,w),(q0,F0,Φ0,b0)) = (
      ∫( (dut⋅ (metric_cf⋅v))*meas_cf )dΩ
    + ∫( (dpt*r)*meas_cf )dΩ
    + ∫( (dBt*w)*meas_cf )dΩ
  )

  c = 0.5 # for jacobian
  jac_x(t,((u,p,B),(q,F,Φ,b)),(du,dp,dB),(v,r,w),(q0,F0,Φ0,b0)) = (
       ∫( du⋅(metric_cf⋅v)*meas_cf )dΩ
     + ∫( (dp*r)*meas_cf )dΩ
     + ∫( (dB*w)*meas_cf )dΩ
    )
  # function jac_prog(dΩ,c)
  #   _jac_prog((t,dt),(u0,h0,B0),(u,h,B),(du,dh,dB),(v,w,r),(b),(F,Φ,q,ω),b3,b1) = (
  #       ∫( du⋅v  )dΩ
  #     + ∫( (c*dt)*(ω*(vecPerp∘(du)⋅v) )  )dΩ
  #     - ∫( ((c*dt)*(1/2))*dB*(∇⋅v) )dΩ
  #     - ∫( ((c*dt)*(1/2)*b1*dh)*(∇⋅v )  )dΩ
  #     + ∫( dh*w   )dΩ
  #     + ∫( (c*dt)*h0*(∇⋅du)*w  )dΩ
  #     + ∫( dB*r )dΩ
  #     + ∫( ((c*dt)*b1*h0)*(∇⋅du)*r )dΩ
  #   )
  # end


  opT = TransientSemilinearFEOperator(mass,res_x,(jac_x,jac_xt),X_prog,Y_prog,assembler=assem_prog)
  opFE = FEOperator(res_y,jac_y,X_diag,Y_diag,assem_diag)
  opDAE = DAEFEOperator(opT,opFE,ls_diag)

  # transient parameters
  _dt = dx(nc(panel_model))*CFL/(p_fe*sqrt(gravity*_H_0))
  Ddt = floor(_dt, sigdigits=1)

  nsteps = _tF/ Ddt
  dt = _tF/floor(nsteps)

  i_am_main(ranks) && println("nsteps = $nsteps")
  i_am_main(ranks) && println("dt = $dt, other dt = $Ddt")

  # solve with SSP RK 3
  ode_solver = RungeKutta(ls_ode,ls_ode,dt,:EXRK_SSP_3_3)
  solT  = solve(ode_solver,opDAE,t0,_tF,xh0)

  ## iterate solution
  it = iterate(solT)

  unwrap_tsw(it,ranks,solT,dir,_tF)


end

function unwrap_tsw(it,ranks,solT,dir,tF,freq=25)
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

  das = FullyAssembledRows()

  # get the ranks to help with storing/saving solution
  ranks = get_ranks(panel_model)

  sim_dir = dir*"/sim_data"
  prog_dir = sim_dir*"/prognostics"
  diag_dir = sim_dir*"/diagnostics"

  vtk_dir = dir*"/vtk_data"
  (i_am_main(ranks) && !isdir(vtk_dir) ) && mkdir(vtk_dir)

  dir_casimirs = dir*"/casimirs"
  (i_am_main(ranks) && !isdir(dir_casimirs)) && mkdir(dir_casimirs)

  lvl = nref(nc(panel_model))

  ## finite element solver
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(das,panel_model)
  dΩ = Measure(Ω_panel,2*(p_fe+1))

  R = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1)
  H = TransientTrialFESpace(R)

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TransientTrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TransientTrialFESpace(V)

  X_prog = TransientMultiFieldFESpace([U,P,P]) # u, p, B
  Y_prog = MultiFieldFESpace([V,Q,Q]) # u, p, B

  X_diag = TransientMultiFieldFESpace([H,U,P,P]) # q, F, Φ, b
  Y_diag = MultiFieldFESpace([R,V,Q,Q]) # q, F, Φ, b

  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  cor_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  gravity = _g

  _Ω_panel = Triangulation(panel_model)
  cell_geo_map = geo_map_func(_Ω_panel)

  labels = ["uh","ph","Bh","bh","qh","Fh","Phih","vort"]
  function make_vtk(t::Float64,xh,yh,cell_geo_map)
    uh,ph,Bh = xh
    qh,Fh,Φh,bh = yh
    vort = qh*ph - cor_cf
    panel_cfs = [covarient_basis_cf⋅uh, ph, Bh, bh, qh, Fh, Φh, vort]

    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(_Ω_panel,vtk_dir*"/solT_$t" * ".vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  function casimirs(xh,yh,dΩ)
    uh,ph,Bh = xh
    qh,Fh,Φh,bh = yh
    vort = qh*ph - cor_cf

    ens = sum(∫( 0.5*(bh*bh*xh[2])*meas_cf  )dΩ)
    energy = sum(∫( (0.5*xh[2]*( xh[1] ⋅(metric_cf⋅xh[1])) + 0.5*xh[2]*xh[3] )*meas_cf )dΩ)
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
  Entropys = Vector{Float64}(undef,length(folders))
  Vorts = Vector{Float64}(undef,length(folders))

  for (i,(f,g)) in enumerate(zip(folders,dfolders))
    t = parse(Float64,f[length(simName)+2:length(f)])

    x =  pload(joinpath(prog_dir,f),ranks)
    xh = FEFunction(X_prog,x)

    y =  pload(joinpath(diag_dir,f),ranks)
    yh = FEFunction(X_diag,y)

    i_am_main(ranks) && println("t = ", t)

    ts[i] = t
    Masss[i], Energys[i], Entropys[i], Vorts[i] = casimirs(xh,yh,dΩ)

    return_vtk && make_vtk(t,xh,yh,cell_geo_map)

    if mod(i,10) == 0
      dxx = dx(panel_model)
      output = @strdict ts Masss Energys Entropys Vorts dxx
      i_am_main(ranks) && safesave(datadir(dir_casimirs, ("casimirs.jld2")), output)
    end

  end

  dxx = dx(panel_model)
  output = @strdict ts Masss Energys Entropys Vorts dxx
  i_am_main(ranks) && safesave(datadir(dir_casimirs, ("casimirs.jld2")), output)

  _make_pvd_distributed(vtk_dir,"solT",1)

end

function convergence_post_process(panel_model,p_fe::Int,dir::String)

  das = FullyAssembledRows()

  # get the ranks to help with storing/saving solution
  ranks = get_ranks(panel_model)

  initial_dir = dir*"/initial_solution"
  final_dir = dir*"/final_solution"

  lvl = nref(nc(panel_model))

  ## finite element solver
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(das,panel_model)
  dΩ_error = Measure(Ω_panel,6*(p_fe+1))

  R = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1)
  H = TransientTrialFESpace(R)

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TransientTrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
  U = TransientTrialFESpace(V)

  X_prog = TransientMultiFieldFESpace([U,P,P]) # u, p, B
  Y_prog = MultiFieldFESpace([V,Q,Q]) # u, p, B

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
  uh,ph,Bh = xh

  x0 =  pload(joinpath(initial_dir,g),ranks)
  xh0 = FEFunction(X_prog,x0)
  uh0,ph0,Bh0 = xh0

  uh_proj = covarient_basis_cf ⋅ uh
  uh0_proj = covarient_basis_cf ⋅ uh0
  e_u = l2( (uh0_proj - uh_proj),meas_cf,dΩ_error)
  e_p = l2((ph0 - ph),meas_cf,dΩ_error)
  e_B = l2((Bh0 - Bh),meas_cf,dΩ_error)

  ### convergence output for DrWatson
  dir_convergence = dir*"/convergence"
  (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

  n = nc(panel_model)
  dxx = dx(panel_model)
  output = @strdict e_u e_p e_B n dxx p_fe lvl t
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("tsw_nref$(lvl)_p$p_fe.jld2")), output)


end
################################################################################
#### Main run for transient solution
################################################################################
function main_transient(distribute,nprocs;
  restart=false,options="",n_ref_lvls=4,p_fe=1,CFL=0.1,ζ=0.0,return_vtk=true,_ε=1e-4,_soft=false)

  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("transient_tsw_equation")

  h = panel_to_cartesian(h₀(ζ))
  vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
  f = panel_to_cartesian(f₀(ζ))
  B = panel_to_cartesian(B₀(ζ))

  ls_diag = CGSolver(JacobiLinearSolver();rtol=1-16,atol=1e-16,verbose=i_am_main(ranks),name="diagnostic_solver")
  ls_ode = CGSolver(JacobiLinearSolver();rtol=1-16,atol=1e-16,verbose=i_am_main(ranks),name="ode_solver")
  lss = (ls_ode,ls_diag)

  omodel = ParametricOctreeDistributedDiscreteModel(ranks; num_initial_uniform_refinements=n_ref_lvls)
  panel_model = omodel.parametric_dmodel

  _dir = datadir("TransientThermalShallowWater_checkpointing")
  (i_am_main(ranks) && !isdir(_dir)) && mkdir(_dir)

  dir = _dir*"/sol_p$(p_fe)_nref$n_ref_lvls"
  (i_am_main(ranks) && !isdir(dir) && return_vtk) && mkdir(dir)

  transient_tsw_solver(panel_model,p_fe,dir,h,vX,f,B,_ε,_soft,CFL,lss,restart)
  convergence_post_process(panel_model,p_fe,dir)
  post_process(panel_model,p_fe,dir,f,return_vtk)

  i_am_main(ranks) && println("--DONE--")
  @test true
end


# MPI.Init()
# nprocs = prod(MPI.Comm_size(MPI.COMM_WORLD))
# ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

# with_mpi() do distribute
#   main_transient(distribute,nprocs;restart=true,n_ref_lvls=3)
# end
