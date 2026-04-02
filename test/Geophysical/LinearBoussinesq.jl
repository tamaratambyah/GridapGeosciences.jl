"""
solve the linear Boussineq equations in 3D in steady form using manufactured solutions
u + f(̂n×u) + ∇ᵧ(φ) - bn̂ = f₁
φ + c² ∇ᵧ⋅u = f₂
b + N² u⋅̂n = f₃
"""

module LinearBoussinesq

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

import GridapGeosciences.Helpers: RADIUS, THICKNESS
THICKNESS
RADIUS

# MPI.Init()
# ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

# dir = datadir("DistributedLinearisedBoussinesq_nprocs$(length(ranks))")
# (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

include("../convergence_tools.jl")
include("../missing_overloads.jl")

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

LH = a_e # m
LV = ztop/THICKNESS
τ = 1/Ωr # s

_d = d/LV
_Lz = Lz/LV
_R = R/LH
_u_0 = u_0*τ/LH
_Ωr = Ωr*τ
_c = c*τ/LH
_N = N*τ
_ztop = ztop/LV


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

function bn(xyz)
  b = b0(xyz)
  n = normal_vec(xyz)
  b*n
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

function un(xyz)
  u = u0(xyz)
  n = normal_vec(xyz)
  u⋅n
end


function omega(xyz)
  x,y,z = xyz
  θϕr   = xyz2θϕr(xyz)
  θ,ϕ,r = θϕr
  2*_Ωr*sin(ϕ)
end


function linear_boussineseq(panel_model::GridapDistributed.GenericDistributedDiscreteModel{3,3},
  p_fe::Int,dir::String,
  h::Function,vX::Function,f::Function,b::Function,_bn::Function,_un::Function,
  ls=LUSolver(),return_vtk=false)

  das =  FullyAssembledRows() # must have FullyAssembledRows to use b_cf on boundary

  ranks = get_ranks(panel_model)

  i_am_main(ranks) && println("Assembly strategy: $das")

  lvl_h = nref(nc_horizontal(panel_model))
  lvl_v = nref(nc_vertical(panel_model))
  i_am_main(ranks) && println("nref_h = $lvl_h; nref_v = $lvl_v; p_fe = $p_fe")

  panel_ids = get_panel_ids(panel_model)
  Ω_panel = Triangulation(das,panel_model)
  dΩ = Measure(Ω_panel,4*(p_fe+1))

  Ω_error = Triangulation(panel_model)
  dΩ_error = Measure(Ω_error,6*p_fe+1)

  covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
  pinvJ_cf = panelwise_cellfield(forward_pinv_jacobian,Ω_panel,panel_ids)

  h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
  u_proj_cf = panelwise_cellfield(projection_v(vX),Ω_panel,panel_ids)
  omega_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
  b_cf = panelwise_cellfield(b,Ω_panel,panel_ids)

  tags = ["bottom_boundary",  "top_boundary"]

  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv,dirichlet_tags=tags)
  U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))

  W = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2,dirichlet_tags=tags)
  B = TrialFESpace(W,b_cf)

  Y = MultiFieldFESpace([V, Q, W])
  X = MultiFieldFESpace([U, P, B])

  u_perp_contra = panelwise_cellfield(contra_v_perp3D(vX),Ω_panel,panel_ids)
  u_perp = covarient_basis_cf ⋅ u_perp_contra

  sgrad_cf = panelwise_cellfield(sgrad(h),Ω_panel,panel_ids)
  sdiv_cf =  panelwise_cellfield(surfdiv(contra_v(vX)),Ω_panel,panel_ids)
  bn_cf = panelwise_cellfield(_bn,Ω_panel,panel_ids)
  un_cf = panelwise_cellfield(_un,Ω_panel,panel_ids)
  g_star_cf = panelwise_cellfield(g_star,Ω_panel,panel_ids)

  # manufacture rhs functions
  rhs_bouyancy = b_cf + _N^2*un_cf
  rhs_pressure = h_cf + _c^2*sdiv_cf
  rhs_vector = u_proj_cf + omega_cf*u_perp + sgrad_cf -bn_cf
  rhs_con_vector = pinvJ_cf ⋅ rhs_vector # exact contravariant component

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

  biform1((u,p,b),(v,q,r)) = ( ∫( (u⋅ (metric_cf⋅v))*meas_cf )dΩ + ∫( ( omega_cf*( (Rperp_cf⋅ u)⋅v))*detg_cf )dΩ
                             - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
                             - ∫( b*(g_star_cf⋅v )*meas_cf )dΩ )
  liform1((v,q,r)) = ∫( rhs_con_vector⋅(metric_cf⋅v)*meas_cf )dΩ

  #### Pressure
  biform2((u,p,b),(v,q,r)) = ∫( (p*q)*meas_cf )dΩ + ∫( _c^2*( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) ) ) )dΩ
  liform2((v,q,r)) =  ∫( (rhs_pressure*q)*meas_cf )dΩ

  #### Bouyancy
  n_cf = CellField(VectorValue(1,0,0),Ω_panel)
  biform3((u,p,b),(v,q,r)) = ∫( (b*r)*meas_cf )dΩ + ∫( _N^2*( r*(g_star_cf⋅u)*meas_cf)   )dΩ
                                                  # ∫( _N^2*( r*( n_cf⋅(metric_cf⋅u))*meas_cf)   )dΩ
  liform3((v,q,r)) =  ∫( (rhs_bouyancy*r)*meas_cf )dΩ


  #### Multifield problem
  assem = SparseMatrixAssembler(X,Y,das)
  biformX((u,p,b),(v,q,r)) = biform1((u,p,b),(v,q,r)) + biform2((u,p,b),(v,q,r)) + biform3((u,p,b),(v,q,r))
  liformX((v,q,r)) = liform1((v,q,r)) + liform2((v,q,r)) + liform3((v,q,r))

  op = AffineFEOperator(biformX,liformX,X,Y,assem)
  A = get_matrix(op)
  b_vec = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b_vec)
  xh = FEFunction(X,x)
  uh,ph,bh = xh

  uh_proj = covarient_basis_cf ⋅ uh
  e_u = l2( (u_proj_cf - uh_proj),meas_cf,dΩ_error) # error in physical velocity u
  e_p = l2((h_cf - ph),meas_cf,dΩ_error) # error in depth
  e_b = l2((b_cf - bh),meas_cf,dΩ_error) # error in bouyancy

  ### solve bouyancy as single field
  u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_contra_h = interpolate(u_contra_cf,U)
  biformB(b,r) = ∫( (b*r)*meas_cf )dΩ
  liformB(r) =  ∫( (rhs_bouyancy*r)*meas_cf )dΩ - ∫( _N^2*( r*( n_cf⋅(metric_cf⋅u_contra_h))*meas_cf)   )dΩ
  op = AffineFEOperator(biformB,liformB,B,W,SparseMatrixAssembler(B,W,das))
  bh = solve(ls,op)
  e_b = l2((b_cf - bh),meas_cf,dΩ_error)


  if return_vtk
    cell_geo_map = geo_map_func(get_panel_ids(Ω_panel))
    panel_cfs = [h_cf, u_proj_cf, b_cf, ph, uh_proj, bh, h_cf-ph, u_proj_cf-uh_proj , b_cf-bh]
    labels = ["p","u_proj", "b", "ph", "uh_proj", "bh", "ep","eu", "eb"]
    cellfields = map((x,y) -> x=>y, labels,panel_cfs)
    writevtk(Ω_panel,dir*"/ambient_model_nrefh$(lvl_h)_nrefv$(lvl_v)_p$p_fe",cellfields=cellfields,append=false,geo_map=cell_geo_map)
  end

  ### convergence output for DrWatson
  dir_convergence = dir*"/convergence"
  (i_am_main(ranks) && !isdir(dir_convergence)) && mkdir(dir_convergence)

  n = nc(panel_model)
  n_h = nc_horizontal(panel_model)
  n_v = _nc_vertical(panel_model)
  dxx = dx(panel_model)
  dxx_horizontal = dx_horizontal(panel_model)
  dxx_vertical = dx_vertical(panel_model)
  output = @strdict e_u e_p e_b n n_h n_v dxx dxx_horizontal dxx_vertical p_fe lvl_h lvl_v
  i_am_main(ranks) && safesave(datadir(dir_convergence, ("linear_boussineseq_nrefh$(lvl_h)_nrefv$(lvl_v)_p$p_fe.jld2")), output)

  return e_u,e_p,e_b
end

################################################################################
#### Auto convergence test
################################################################################
function main(distribute,nprocs)
  ranks = distribute(LinearIndices((nprocs,)))

  i_am_main(ranks) && println("--START--")
  i_am_main(ranks) && println("Auto conference test: Linear Boussineq")

  dir = datadir("LinearBoussineseqConvergence")
  (i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

  n_ref_lvls = 3
  ps = [1]
  ls = LUSolver()

  models = get_3D_octree_refined_models(ranks,n_ref_lvls)

  h = panel_to_cartesian(p0)
  vX = panel_to_cartesian(tangent_vec(u0))
  f = panel_to_cartesian(omega)
  b = panel_to_cartesian(b0)
  _bn = panel_to_cartesian(bn)
  _un = panel_to_cartesian(un)

  p_convergence_test(ranks,ps,models,linear_boussineseq,dir,h,vX,f,b,_bn,_un,ls,true)


  i_am_main(ranks) && println("--DONE--")
end



end ## module
