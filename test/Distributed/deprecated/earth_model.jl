using MPI
using PartitionedArrays
using GridapGeosciences
using GridapGeosciences.Distributed
using GridapP4est
using Gridap
using GridapDistributed

using DrWatson

using Gridap.Algebra
include("../helpers.jl")
include("../missing_overloads.jl")

MPI.Init()
ranks = distribute_with_mpi(LinearIndices((prod(MPI.Comm_size(MPI.COMM_WORLD)),)))

dir = datadir("EarthModel")
(i_am_main(ranks) && !isdir(dir)) && mkdir(dir)

import GridapGeosciences.Helpers: RADIUS, THICKNESS
THICKNESS
RADIUS

num_horizontal_uniform_refinements = 3
num_vertical_uniform_refinements = 3
o3model = GridapGeosciences.Distributed.Parametric3DOctreeDistributedDiscreteModel(ranks;
                                           num_horizontal_uniform_refinements=num_horizontal_uniform_refinements,
                                           num_vertical_uniform_refinements=num_vertical_uniform_refinements);

panel_model = o3model.parametric_dmodel
Ω_panel =  Triangulation(panel_model)
panel_ids = get_panel_ids(Ω_panel)
panel_ids = get_panel_ids(panel_model)


gibson_dx = sqrt(4*π*(6.36e6/125)^2/24576)
gibson_dz = 10e3/64
gibson_dx/gibson_dz

ratio = 7.5#7.37

_gibson_dx = sqrt(4*π*1^2/24576)
_gibson_tt = _gibson_dx/ratio* 64
_gibson_dz = _gibson_tt/64
_gibson_dx/_gibson_dz

horizontal = 4*π*RADIUS^2/(nc_horizontal(panel_model)*6)
dxx =  sqrt(horizontal)
tt = dxx/ratio*_nc_vertical(panel_model)
tt = THICKNESS
dzz =tt/_nc_vertical(panel_model)
dxx/dzz



## Try plotting in distributed
cell_geo_map = geo_map_func(Ω_panel)
writevtk(Ω_panel,dir*"/model_radius",append=false,geo_map=cell_geo_map)


######### linear Boussineq
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
TF = (3600*24) # s

LH = a_e # m
LV = ztop/THICKNESS #m
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

h = panel_to_cartesian(p0)
vX = panel_to_cartesian(tangent_vec(u0))
f = panel_to_cartesian(omega)
b = panel_to_cartesian(b0)
_bn = panel_to_cartesian(bn)
_un = panel_to_cartesian(un)

p_fe = 1

das =  FullyAssembledRows()

panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(das,panel_model)
dΩ = Measure(Ω_panel,4*(p_fe+1))

Ω_error = Triangulation(panel_model)
dΩ_error = Measure(Ω_error,6*p_fe+1)

covariant_basis_cf = panelwise_cellfield(covariant_basis,Ω_panel,panel_ids)
pinvJ_cf = panelwise_cellfield(forward_pinv_jacobian,Ω_panel,panel_ids)

h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
u_proj_cf = panelwise_cellfield(projection_v(vX),Ω_panel,panel_ids)
omega_cf = panelwise_cellfield(f,Ω_panel,panel_ids)
b_cf = panelwise_cellfield(b,Ω_panel,panel_ids)



cell_geo_map = geo_map_func(get_panel_ids(Ω_panel))
panel_cfs = [h_cf, u_proj_cf, b_cf,omega_cf ]
labels = ["p","u_proj", "b",  "omega"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/ambient_model_nrefh",cellfields=cellfields,append=false,geo_map=cell_geo_map)


## solve
tags = ["bottom_boundary",  "top_boundary"]



  Q = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
  P = TrialFESpace(Q)

  V = TestFESpace(Ω_panel, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv,dirichlet_tags=tags)
  U = TrialFESpace(V,VectorValue(0.0,0.0,0.0))


  W = TestFESpace(Ω_panel, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2,dirichlet_tags=tags)
  B = TrialFESpace(W,b_cf)

  Y = MultiFieldFESpace([V, Q, W])
  X = MultiFieldFESpace([U, P, B])

  u_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  u_perp_contra = panelwise_cellfield(contra_v_perp3D(vX),Ω_panel,panel_ids)
  u_perp = covariant_basis_cf ⋅ u_perp_contra

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
  # biform3((u,p,b),(v,q,r)) = ∫( (b*r)*meas_cf )dΩ + ∫( _N^2*( r*(g_star_cf⋅u)*meas_cf)   )dΩ
  biform3((u,p,b),(v,q,r)) = ∫( (b*r)*meas_cf )dΩ + ∫( _N^2*( r*( n_cf⋅(metric_cf⋅u))*meas_cf)   )dΩ
  liform3((v,q,r)) =  ∫( (rhs_bouyancy*r)*meas_cf )dΩ


  #### Multifield problem
  assem = SparseMatrixAssembler(X,Y,das)
  biformX((u,p,b),(v,q,r)) = biform1((u,p,b),(v,q,r)) + biform2((u,p,b),(v,q,r)) + biform3((u,p,b),(v,q,r))
  liformX((v,q,r)) = liform1((v,q,r)) + liform2((v,q,r)) + liform3((v,q,r))



  ls = LUSolver()
  op = AffineFEOperator(biformX,liformX,X,Y,assem)
  A = get_matrix(op)
  b_vec = get_vector(op)
  ns = numerical_setup(symbolic_setup(ls,A),A)
  x = allocate_in_domain(A); fill!(x,0.0)
  solve!(x,ns,b_vec)
  xh = FEFunction(X,x)
  uh,ph,bh = xh

  uh_proj = covariant_basis_cf ⋅ uh
  e_u = l2( (u_proj_cf - uh_proj),meas_cf,dΩ_error) # error in physical velocity u
  # e_u = l2( (u_contra_cf - uh),dΩ_error) # error in contr u
  e_p = l2((h_cf - ph),meas_cf,dΩ_error) # error in depth
  e_b = l2((b_cf - bh),meas_cf,dΩ_error) # error in bouyancy

  u_contra_h = interpolate(u_contra_cf,U)
  #### Pressure
  biformP(p,q) = ∫( (p*q)*meas_cf )dΩ
  liformP(q) =  ∫( (rhs_pressure*q)*meas_cf )dΩ - ∫( _c^2*( q*(uh⋅grad_meas_cf + meas_cf*(∇⋅uh) ) ) )dΩ
  op = AffineFEOperator(biformP,liformP,P,Q,SparseMatrixAssembler(P,Q,das))
  ph = solve(ls,op)
  e_p = l2((h_cf - ph),meas_cf,dΩ_error)

  ### bouyancy
  n_cf = CellField(VectorValue(1,0,0),Ω_panel)
  # biform3((u,p,b),(v,q,r)) = ∫( (b*r)*meas_cf )dΩ + ∫( _N^2*( r*(g_star_cf⋅u)*meas_cf)   )dΩ
  biformB(b,r) = ∫( (b*r)*meas_cf )dΩ
  liformB(r) =  ∫( (rhs_bouyancy*r)*meas_cf )dΩ - ∫( _N^2*( r*( n_cf⋅(metric_cf⋅u_contra_h))*meas_cf)   )dΩ
  op = AffineFEOperator(biformB,liformB,B,W,SparseMatrixAssembler(B,W,das))
  bh = solve(ls,op)
  e_b = l2((b_cf - bh),meas_cf,dΩ_error)

 cell_geo_map = geo_map_func(get_panel_ids(Ω_panel))
panel_cfs = [h_cf, u_proj_cf, b_cf,omega_cf, meas_cf, uh_proj,ph,bh, u_contra_cf - uh,h_cf - ph,b_cf - bh]
labels = ["p","u_proj", "b",  "omega", "g", "uh_proj", "ph","bh", "eu", "ep", "eb"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/ambient_model_nrefh",cellfields=cellfields,append=false,geo_map=cell_geo_map)
