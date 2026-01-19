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

include("../convergence_tools.jl")
include("Williamson5Test.jl")


function transient_shallow_water_solver_3D(
  panel_model::GridapDistributed.DistributedDiscreteModel{3,3},
  p_fe::Int,_dir::String,h::Function,vX::Function,f::Function,b::Function,
  lss=(LUSolver(),LUSolver()),CFL=0.1,return_vtk=false)

  ls_ode, ls_diag = lss

  # get the ranks to help with storing/saving solution
  ranks = get_ranks(panel_model)

  lvl_h = nref(nc_horizontal(panel_model))
  lvl_v = nref(nc_vertical(panel_model))
  i_am_main(ranks) && println("nref_h = $lvl_h; nref_v = $lvl_v; p_fe = $p_fe")

  lvl = lvl_h
  i_am_main(ranks) && println("nlevl = $lvl")

  dir = _dir*"/sol_p$(p_fe)_nref$lvl"
  (i_am_main(ranks) && !isdir(dir) && return_vtk) && mkdir(dir)

  dir_latlon = _dir*"/latlon_sol_p$(p_fe)_nref$lvl"
  (i_am_main(ranks) && !isdir(dir_latlon) && return_vtk) && mkdir(dir_latlon)


  ## finite element solver
  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(panel_model)
  dΩ = Measure(Ω_panel,4*(p_fe+1))


  tags = ["top_boundary", "bottom_boundary"]

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv,dirichlet_tags=tags)
  U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))

  R = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1)
  H = TransientTrialFESpace(R)

  Y_prog = MultiFieldFESpace([V, Q]) # u, p
  X_prog = TransientMultiFieldFESpace([U, P]) # u, p

  X_diag = TransientMultiFieldFESpace([H,U,P]) # q, F, Φ
  Y_diag = MultiFieldFESpace([R,V,Q]) # q, F, Φ

  # initial conditions
  h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
  b_cf = panelwise_cellfield(b,Ω_panel,panel_ids)
  h_h = interpolate(h_cf-b_cf,P)
  u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_contra_h = interpolate(u_contra_cf,U)
  xh0 = interpolate([u_contra_h,h_h],X_prog)

  cor_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  gravity = _g

  # weak forms
  detg_cf = panelwise_cellfield(detg,Ω_panel,panel_ids)
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)
  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

  #### DIAGNOSTIC VARIABLES

  # vorticity
  function perp_matrix_3D(p::Int,αβ)
    @check length(αβ) == 3
    _m = metric(p,αβ)
    M = [ _m[2,2] _m[2,3] ; _m[3,2] _m[3,3] ]
    m = TensorValue(M)
    TensorValue{2,2}( -m[1,2], m[1,1], -m[2,2], m[1,2] )
  end
  perp_matrix_3D(p) = αβ -> perp_matrix_3D(p,αβ)

  function extract_2D(u, m)
    u2D = VectorValue(u[2],u[3])
    v = m⋅u2D
    return VectorValue(0,v[1],v[2])
  end

  perp_matrix_cf = panelwise_cellfield(perp_matrix_3D,Ω_panel,panel_ids)
  resq(((u,p),(q,F,Φ)),(w,v,ψ)) = ∫( q*p*w*meas_cf  )dΩ - ∫( cor_cf*w*meas_cf  )dΩ - ∫( (( extract_2D∘(u,perp_matrix_cf)) )⋅∇(w)  )dΩ

  # mass flux
  resF(((u,p),(q,F,Φ)),(w,v,ψ)) = ∫( (F⋅ (metric_cf⋅v))*meas_cf )dΩ - ∫( p*(u⋅(metric_cf⋅v))*meas_cf   )dΩ

  # Bernoulli potential
  resΦ(((u,p),(q,F,Φ)),(w,v,ψ)) = ∫( Φ*ψ*meas_cf  )dΩ - ∫( gravity*(p+b_cf)*ψ*meas_cf  )dΩ - ∫( 0.5*( u ⋅(metric_cf⋅u) )ψ*meas_cf  )dΩ

  res_y(t,((u,p),(q,F,Φ)),(w,v,ψ)) = resq(((u,p),(q,F,Φ)),(w,v,ψ)) + resF(((u,p),(q,F,Φ)),(w,v,ψ)) + resΦ(((u,p),(q,F,Φ)),(w,v,ψ))
  jac_y(t,((u,p),(q,F,Φ)),(dq,dF,dΦ),(w,v,ψ)) = ∫( dq*p*w*meas_cf  )dΩ + ∫( (dF⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( dΦ*ψ*meas_cf  )dΩ

  _res_y((q,F,Φ),(w,v,ψ))  = res_y(0.0,(xh0,(q,F,Φ)),(w,v,ψ))
  _jac_y((q,F,Φ),(dq,dF,dΦ),(w,v,ψ)) = jac_y(0.0,(xh0,(q,F,Φ)),(dq,dF,dΦ),(w,v,ψ))
  _opFE = FEOperator(_res_y,_jac_y,X_diag,Y_diag)
  nls = GridapSolvers.NonlinearSolvers.NewtonSolver(ls_diag,verbose=i_am_main(ranks))
  qh,Fh,Φh = solve(nls,_opFE)
  vort = qh*xh0[2] - cor_cf

  #### PROGNOSTIC VARIABLES

  # equation for depth and velocity:
  mass(t,(dut,dpt),(v,r)) = ∫( (dut⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( (dpt*r)*meas_cf )dΩ
  res_p(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) = ∫( r*(F⋅grad_meas_cf + meas_cf*(∇⋅F) )  )dΩ

  #### Velocity
  Aperp = [0 0 0
            0 0 -1
            0 1 0]
  Rperp = TensorValue(Aperp)
  Rperp_cf = CellField(Rperp,Ω_panel)

  res_u(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) = (
                                ∫( q*( (Rperp_cf⋅ F)⋅v)*detg_cf  )dΩ
                              + ∫( -τ*( (q-q0)/dt )*( (Rperp_cf⋅ F)⋅v)*detg_cf  )dΩ
                              + ∫( -τ*(u⋅∇(q))*( (Rperp_cf⋅ F)⋅v)*detg_cf  )dΩ
                              - ∫( Φ*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
                    )

  res_x(t,((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) = res_u(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0)) + res_p(((u,p),(q,F,Φ)),(v,r),(q0,F0,Φ0))
  jac_x(t,((u,p),(q,F,Φ)),(du,dp),(v,r),(q0,F0,Φ0)) = ∫( -τ*(du⋅∇(q))*( (Rperp_cf⋅ F)⋅v)*detg_cf  )dΩ
  jac_xt(t,((u,p),(q,F,Φ)),(dut,dpt),(v,r),(q0,F0,Φ0)) =  ∫( (dut⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( (dpt*r)*meas_cf )dΩ


  opT = TransientSemilinearFEOperator(mass,res_x,(jac_x,jac_xt),X_prog,Y_prog)
  opFE = FEOperator(res_y,jac_y,X_diag,Y_diag)
  opDAE = DAEFEOperator(opT,opFE,ls_diag)

  t0, tF = 0.0, _tF
  # _dt = dx(panel_model)*CFL/(p_fe*sqrt(gravity*_H_0))
  _dt = dx(nc(panel_model))*CFL/(p_fe*sqrt(gravity*_H_0))
  dt = floor(_dt, sigdigits=1)

  τ = dt/2

  ode_solver = RungeKutta(ls_ode,ls_ode,dt,:EXRK_SSP_3_3)

  solT  = solve(ode_solver,opDAE,t0,tF,xh0)
  it = iterate(solT)



  cell_geo_map = geo_map_func(Ω_panel)
  latlon_cell_geo_map = latlon_geo_map_func(Ω_panel)
  owned_panel_ids = get_owned_panel_ids(panel_model)
  if return_vtk
    panel_cfs = [covarient_basis_cf⋅xh0[1], xh0[2],qh,Fh,Φh,vort,b_cf, owned_panel_ids]
    cellfields = map((x,y) -> x=>y, ["uh","ph","qh","Fh","Phih","vort","bt","pid"],panel_cfs)
    writevtk(Ω_panel,dir*"/solT_0.vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map,order=2)
    writevtk(Ω_panel,dir_latlon*"/latlon_solT_0.vtu", cellfields=cellfields,append=false,geo_map=latlon_cell_geo_map,order=2)
  end


  counter = 1
  while !isnothing(it)
    data, state = it
    t, xh = data
    odeopcache = state[2][5][2]
    yh = odeopcache.diagnostics

    uh,ph = xh
    qh,Fh,Φh = yh

    vort = qh*ph - cor_cf
    i_am_main(ranks) && println(t)

    if return_vtk
      panel_cfs = [covarient_basis_cf⋅uh, ph,qh,Fh,Φh,vort,owned_panel_ids]
      cellfields = map((x,y) -> x=>y, ["uh","ph","qh","Fh","Phih","vort","pid"],panel_cfs)
      writevtk(Ω_panel,dir*"/solT_$t.vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map,order=2)
      writevtk(Ω_panel,dir_latlon*"/latlon_solT_$t.vtu", cellfields=cellfields,append=false,geo_map=latlon_cell_geo_map,order=2)

    end
    counter = counter + 1
    it = iterate(solT, state)
  end

  # _make_pvd_distributed(dir,"solT",1)

end


################################################################################
#### Main run for transient solution
################################################################################
function main_transient(distribute,nprocs;n_ref_lvls=4,p_fe=1,CFL=0.1,ζ=0.0,return_vtk=true)
  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("transient_shallow_water_3D")

  h = panel_to_cartesian(h₀(ζ))
  vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
  f = panel_to_cartesian(f₀(ζ))
  b = panel_to_cartesian(topography)

  o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
  num_horizontal_uniform_refinements=n_ref_lvls,num_vertical_uniform_refinements=0)
  panel_model = o3model.parametric_dmodel

  ls_diag = CGSolver(JacobiLinearSolver();rtol=1-12,verbose=i_am_main(ranks),name="diagnostic_solver")
  ls_ode = CGSolver(JacobiLinearSolver();rtol=1-12,verbose=i_am_main(ranks),name="ode_solver")

  lss = (ls_ode,ls_diag)

  dir = datadir("TransientShallowWater_W5_3D")
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  transient_shallow_water_solver_3D(panel_model,p_fe,dir,h,vX,f,b,lss,CFL,return_vtk)


  i_am_main(ranks) && println("--DONE--")

end


# with_mpi() do distribute
#   main_transient(distribute,2;n_ref_lvls=3)
# end

# dir = datadir("TransientShallowWater_W5_3D/sol_p1_nref3")
# _make_pvd_distributed(dir,"solT",1)
