using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra
using GridapPETSc
using GridapGeosciences
using Test
using CSV
using DataFrames
using MPI

include("../convergence_tools.jl")
include("../output_tools.jl")
include("ThermogeostrophicBalanceTest.jl")

# ζ = 0.0
# h = panel_to_cartesian(h₀(ζ))
# vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
# f = panel_to_cartesian(f₀(ζ))
# B = panel_to_cartesian(B₀(ζ))


# n_ref_lvls = 4
# p_fe = 1
# lss = (LUSolver(),LUSolver())
# return_vtk = true
# CFL = 0.1
# models  = get_refined_models(n_ref_lvls)

# panel_model = models[1]

# _dir = datadir("TransientICTest")

# tsw_solver(panel_model,p_fe,_dir,h,vX,f,B,lss,CFL,return_vtk)

# d = _dir*"/latlon_sol_p1_nref4"
# make_pvd(d,"latlon_solT",1)

function my_mean( Bu_n::SkeletonPair)
  plus  = ( Bu_n.plus)
  minus = ( Bu_n.minus)
  0.5*( plus - minus  )
end

function upwinding_sign(Fn)
  c = 0.0

  if Fn < 0.0
    c = -0.5
  elseif Fn > 0.0
    c = 0.5
  end
  return c
end


function tsw_solver(
  panel_model::Union{<:DiscreteModel{2,2},<:GridapDistributed.DistributedDiscreteModel{2,2}},
  p_fe::Int,_dir::String,h::Function,vX::Function,f::Function,B::Function,
  lss=(LUSolver(),LUSolver()),CFL=0.1,return_vtk=false)

  das = FullyAssembledRows()

  ls_ode, ls_diag = lss

  # get the ranks to help with storing/saving solution
  ranks = get_ranks(panel_model)
  i_am_main(ranks) && println("Assembly strategy: $das")


  lvl = nref(nc(panel_model))
  i_am_main(ranks) && println("nlevl = $lvl")

  dir = _dir*"/sol_p$(p_fe)_nref$lvl"
  (i_am_main(ranks) && !isdir(dir) && return_vtk) && mkdir(dir)

  dir_latlon = _dir*"/latlon_sol_p$(p_fe)_nref$lvl"
  (i_am_main(ranks) && !isdir(dir_latlon) && return_vtk) && mkdir(dir_latlon)


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

  X_diag = TransientMultiFieldFESpace([H,U,P,P,P]) # q, F, Φ, T, b
  Y_diag = MultiFieldFESpace([R,V,Q,Q,Q]) # q, F, Φ, T, b


  ## initial conditions
  covariant_basis_cf = panelwise_cellfield(covariant_basis,Ω_panel,panel_ids)
  u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_contra_h = interpolate(u_contra_cf,U)
  u_proj_h = covariant_basis_cf ⋅ u_contra_h

  h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
  B_cf = panelwise_cellfield(B,Ω_panel,panel_ids)
  h_h = interpolate(h_cf,P)
  B_h = interpolate(B_cf,P)

  xh0 = interpolate_everywhere([u_contra_h,h_h,B_h],X_prog(0.0))

  cor_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  gravity = _g


  # mectrics required in weak forms
  metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
  meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
  grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)
  meas_cf_skel = panelwise_cellfield(sqrtg,Λ)


  #### DIAGNOSTIC VARIABLES
  assem_diag = SparseMatrixAssembler(X_diag,Y_diag,das)

  # vorticity
  perp_matrix_cf = panelwise_cellfield(perp_matrix,Ω_panel,panel_ids)
  resq(((u,p,B),(q,F,Φ,T,b)),(w,v,ψ,s,r)) = ∫( q*p*w*meas_cf  )dΩ - ∫( cor_cf*w*meas_cf  )dΩ - ∫( (perp_matrix_cf⋅u)⋅∇(w)  )dΩ

  # mass flux
  resF(((u,p,B),(q,F,Φ,T,b)),(w,v,ψ,s,r)) = ∫( (F⋅ (metric_cf⋅v))*meas_cf )dΩ - ∫( p*(u⋅(metric_cf⋅v))*meas_cf   )dΩ

  # Bernoulli potential
  resΦ(((u,p,B),(q,F,Φ,T,b)),(w,v,ψ,s,r)) = ∫( Φ*ψ*meas_cf  )dΩ - ∫( 0.5*B*ψ*meas_cf  )dΩ - ∫( 0.5*( u ⋅(metric_cf⋅u) )ψ*meas_cf  )dΩ

  # Temperature
  resT(((u,p,B),(q,F,Φ,T,b)),(w,v,ψ,s,r)) = ∫( T*s*meas_cf  )dΩ - ∫( 0.5*p*s*meas_cf  )dΩ

  # Bouyancy
  resb(((u,p,B),(q,F,Φ,T,b)),(w,v,ψ,s,r)) = ∫( b*p*r*meas_cf  )dΩ - ∫( B*r*meas_cf  )dΩ


  res_y(t,((u,p,B),(q,F,Φ,T,b)),(w,v,ψ,s,r)) = (
      resq(((u,p,B),(q,F,Φ,T,b)),(w,v,ψ,s,r))
    + resF(((u,p,B),(q,F,Φ,T,b)),(w,v,ψ,s,r))
    + resΦ(((u,p,B),(q,F,Φ,T,b)),(w,v,ψ,s,r))
    + resT(((u,p,B),(q,F,Φ,T,b)),(w,v,ψ,s,r))
    + resb(((u,p,B),(q,F,Φ,T,b)),(w,v,ψ,s,r))
  )
  jac_y(t,((u,p,B),(q,F,Φ,T,b)),(dq,dF,dΦ,dT,db),(w,v,ψ,s,r)) = (
      ∫( dq*p*w*meas_cf  )dΩ
    + ∫( (dF⋅ (metric_cf⋅v))*meas_cf )dΩ
    + ∫( dΦ*ψ*meas_cf  )dΩ
    + ∫( dT*s*meas_cf  )dΩ
    + ∫( db*p*r*meas_cf  )dΩ
  )

  _res_y((q,F,Φ,T,b),(w,v,ψ,s,r))  = res_y(0.0,(xh0,(q,F,Φ,T,b)),(w,v,ψ,s,r))
  _jac_y((q,F,Φ,T,b),(dq,dF,dΦ,dT,db),(w,v,ψ,s,r)) = jac_y(0.0,(xh0,(q,F,Φ,T,b)),(dq,dF,dΦ,dT,db),(w,v,ψ,s,r))
  _opFE = FEOperator(_res_y,_jac_y,X_diag,Y_diag,assem_diag)
  nls = GridapSolvers.NonlinearSolvers.NewtonSolver(ls_diag,verbose=i_am_main(ranks))
  qh,Fh,Φh,Th,bh = solve(nls,_opFE)
  vort = qh*xh0[2] - cor_cf




  #### PROGNOSTIC VARIABLES
  assem_prog = SparseMatrixAssembler(X_prog,Y_prog,das)

  # equation for depth and velocity:
  mass(t,(dut,dpt,dBt),(v,r,w)) = (
      ∫( (dut⋅ (metric_cf⋅v))*meas_cf )dΩ
    + ∫( (dpt*r)*meas_cf )dΩ
    + ∫( (dBt*w)*meas_cf )dΩ
  )

  res_p(((u,p,B),(q,F,Φ,T,b)),(v,r,w),(q0,F0,Φ0,T0,b0)) = ∫( r*(F⋅grad_meas_cf + meas_cf*(∇⋅F) )  )dΩ

  res_u(((u,p,B),(q,F,Φ,T,b)),(v,r,w),(q0,F0,Φ0,T0,b0)) = (
            ∫( q*( (perp_matrix_cf⋅F) ⋅(metric_cf ⋅v))   )dΩ
          - ∫( Φ*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
          + ∫( 0.5*(b*(∇(T)⋅v) )*meas_cf )dΩ
          + ∫( -0.5*(b*T)*(v⋅grad_meas_cf + meas_cf*(∇⋅v) )  )dΩ
          + ∫( -0.5*(T*(∇(b)⋅v) )*meas_cf )dΩ
      )

  u_s1(((u,p,B),(q,F,Φ,T,b)),(v,r,w)) = (
      ∫( -0.5*my_mean((v*b)⋅n_Λ)*jump(T)*meas_cf_skel.plus   )dΛ
    + ∫( 0.5*my_mean((v*T)⋅n_Λ)*jump(b)*meas_cf_skel.plus   )dΛ
  )

  u_s2(((u,p,B),(q,F,Φ,T,b)),(v,r,w)) = ∫( -0.5*( (upwinding_sign∘((F⋅ n_Λ).plus))*(v⋅n_Λ).plus )*jump(b)*jump(T)*meas_cf_skel.plus   )dΛ


  res_B(((u,p,B),(q,F,Φ,T,b)),(v,r,w),(q0,F0,Φ0,T0,b0)) = (
      ∫( -0.5*(b*(∇(w)⋅F) )*meas_cf )dΩ
    + ∫( 0.5*(b*w)*(F⋅grad_meas_cf + meas_cf*(∇⋅F) )  )dΩ
    + ∫( 0.5*(w*(∇(b)⋅F) )*meas_cf )dΩ
  )



  B_s1(((u,p,B),(q,F,Φ,T,b)),(v,r,w)) = (
      ∫( 0.5*my_mean((F*b)⋅n_Λ)*jump(w)*meas_cf_skel.plus   )dΛ
    + ∫( 0.5*my_mean((F*w)⋅n_Λ)*jump(b)*meas_cf_skel.plus   )dΛ
  )

  B_s2(((u,p,B),(q,F,Φ,T,b)),(v,r,w)) = ∫( 0.5*( (upwinding_sign∘((F⋅ n_Λ).plus))*(F⋅n_Λ).plus )*jump(b)*jump(w)*meas_cf_skel.plus   )dΛ


  res_x(t,((u,p,B),(q,F,Φ,T,b)),(v,r,w),(q0,F0,Φ0,T0,b0)) = (
      res_u(((u,p,B),(q,F,Φ,T,b)),(v,r,w),(q0,F0,Φ0,T0,b0))
    + u_s1(((u,p,B),(q,F,Φ,T,b)),(v,r,w))
    + u_s2(((u,p,B),(q,F,Φ,T,b)),(v,r,w))
    + res_p(((u,p,B),(q,F,Φ,T,b)),(v,r,w),(q0,F0,Φ0,T0,b0))
    + res_B(((u,p,B),(q,F,Φ,T,b)),(v,r,w),(q0,F0,Φ0,T0,b0))
    + B_s1(((u,p,B),(q,F,Φ,T,b)),(v,r,w))
    + B_s2(((u,p,B),(q,F,Φ,T,b)),(v,r,w))
  )
  jac_x(t,((u,p,B),(q,F,Φ,T,b)),(du,dp,dB),(v,r,w),(q0,F0,Φ0,T0,b0)) =  ∫( VectorValue(0,0)⋅(du⋅v) + 0*dp*r + 0*dB*w   )dΩ
  jac_xt(t,((u,p,B),(q,F,Φ,T,b)),(dut,dpt,dBt),(v,r,w),(q0,F0,Φ0,T0,b0)) = (
      ∫( (dut⋅ (metric_cf⋅v))*meas_cf )dΩ
    + ∫( (dpt*r)*meas_cf )dΩ
    + ∫( (dBt*w)*meas_cf )dΩ
  )




  opT = TransientSemilinearFEOperator(mass,res_x,(jac_x,jac_xt),X_prog,Y_prog,assembler=assem_prog)
  opFE = FEOperator(res_y,jac_y,X_diag,Y_diag,assem_diag)
  opDAE = DAEFEOperator(opT,opFE,ls_diag)

  t0, tF = 0.0, _tF
  # _dt = dx(panel_model)*CFL/(p_fe*sqrt(gravity*_H_0))
  _dt = dx(nc(panel_model))*CFL/(p_fe*sqrt(gravity*_H_0))
  dt = floor(_dt, sigdigits=1)


  ode_solver = RungeKutta(ls_ode,ls_ode,dt,:EXRK_SSP_3_3)

  solT  = solve(ode_solver,opDAE,t0,tF,xh0)
  it = iterate(solT)

  Ω_error = Triangulation(panel_model)
  dΩ_error = Measure(Ω_error,6*p_fe+1)
  Es_u = Float64[]
  Es_p = Float64[]
  Es_B = Float64[]
  e_u,e_p,e_B = 0.0, 0.0, 0.0
  push!(Es_u,e_u)
  push!(Es_p,e_p)
  push!(Es_B,e_B)

  cell_geo_map = geo_map_func(get_panel_ids(Ω_error))
  latlon_cell_geo_map = latlon_geo_map_func(get_panel_ids(Ω_error))
  if return_vtk
    panel_cfs = [covariant_basis_cf⋅xh0[1], xh0[2], xh0[3],  qh,covariant_basis_cf⋅Fh,Φh,vort,Th,bh]
    cellfields = map((x,y) -> x=>y, ["uh","ph", "Bh",  "qh","Fh","Phih","vort", "Th","bh"],panel_cfs)
    writevtk(Ω_panel,dir*"/solT_0.vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)
    writevtk(Ω_panel,dir_latlon*"/latlon_solT_0.vtu", cellfields=cellfields,append=false,geo_map=latlon_cell_geo_map)
  end


  counter = 1
  while !isnothing(it)
    data, state = it
    t, xh = data
    odeopcache = state[2][5][2]
    yh = odeopcache.diagnostics

    uh,ph,Bh = xh
    qh,Fh,Φh,Th,bh = yh

    vort = qh*ph - cor_cf
    i_am_main(ranks) && println(t)

    uh_proj = covariant_basis_cf ⋅ uh
    e_u = l2( (u_proj_h - uh_proj),meas_cf,dΩ_error)
    e_p = l2((h_cf - ph),meas_cf,dΩ_error)
    e_B = l2((B_cf - Bh),meas_cf,dΩ_error)
    push!(Es_u,e_u)
    push!(Es_p,e_p)
    push!(Es_B,e_B)

    if return_vtk  && (mod(counter,50) == 0)
      panel_cfs = [covariant_basis_cf⋅uh, ph, Bh, bh,  qh,Fh,Φh,vort]
      cellfields = map((x,y) -> x=>y, ["uh","ph", "Bh", "bh", "qh","Fh","Phih","vort"],panel_cfs)
      writevtk(Ω_panel,dir*"/solT_$t.vtu", cellfields=cellfields,append=false,geo_map=cell_geo_map)
      writevtk(Ω_panel,dir_latlon*"/latlon_solT_$t.vtu", cellfields=cellfields,append=false,geo_map=latlon_cell_geo_map)
    end

    counter = counter + 1
    it = iterate(solT, state)
  end

  ### convergence output for DrWatson
  dir_convergence = _dir*"/convergence"
  (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

  n = nc(panel_model)
  dxx = dx(panel_model)
  output = @strdict Es_u Es_p Es_B n dxx p_fe lvl
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("shallow_water_nref$(lvl)_p$p_fe.jld2")), output)

  return Es_u, Es_p, Es_B

end



## helper function to return errors
function transient_tsw_errors(panel_model,p_fe::Int,dir::String,
  h::Function,vX::Function,f::Function,B::Function,lss=(LUSolver(),LUSolver()),CFL=0.1,return_vtk=false)
  Es_u,Es_p,Es_B  = tsw_solver(panel_model,p_fe,dir,h,vX,f,B,lss,CFL,return_vtk)
  return minimum(Es_p[end-10:end]),minimum(Es_u[end-10:end]),minimum(Es_B[end-10:end])
end


################################################################################
#### Auto convergence test
################################################################################
function main(distribute,nprocs;octree=false)
  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("Auto conference test: Thermal Shallow Water")

  n_ref_lvls = 5
  ps = [1]
  lss = (LUSolver(),LUSolver())
  return_vtk = true
  ζ = 0.0

  # ls_diag = CGSolver(JacobiLinearSolver();rtol=1-12,verbose=i_am_main(ranks),name="diagnostic_solver")
  # ls_ode = CGSolver(JacobiLinearSolver();rtol=1-8,verbose=i_am_main(ranks),name="ode_solver")
  # lss = (ls_ode,ls_diag)

  CFL = 0.1

  dir = foldername("ThermalShallowWaterConvergence_upwinded",octree,false)
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  models = get_models(ranks,nprocs,n_ref_lvls;threedims=false,octree=octree)


  h = panel_to_cartesian(h₀(ζ))
  vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
  f = panel_to_cartesian(f₀(ζ))
  B = panel_to_cartesian(B₀(ζ))

  p_convergence_test(ranks,ps,models,transient_tsw_errors,dir,h,vX,f,B,lss,CFL,return_vtk)


  i_am_main(ranks) && println("WARNING! Error output is [p,u,B]")
  i_am_main(ranks) && println("--DONE--")

end



# with_mpi() do distribute
#   main(distribute,2;octree=true)
# end
