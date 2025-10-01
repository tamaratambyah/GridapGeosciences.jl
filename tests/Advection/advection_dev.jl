#### Linear advection equation
#### solve with SUPG as per Brooks & Hughes 1982 paper
#### Replicate tests in Lauritzen2012 paper with time dependent velocity

using DrWatson
using Gridap
using GridapDistributed
using PartitionedArrays
using MPIPreferences

using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers
MPIPreferences.use_jll_binary()

using GridapGeosciences
using GridapSolvers
using Plots, LaTeXStrings
using GridapPETSc


include("Lauritzen_functions.jl")
include("../Distributed/distributed_model.jl")

dir = datadir("Transient_advection_test_supg_debug2")
!isdir(dir) && mkdir(dir)

nprocs = 6
ranks  = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end

panel_model = get_distributed_panel_model(ranks,nprocs,6)


p_fe = 1
CFL = 0.1


u0(XYZ) = cosine_bell(XYZ)
u = panel_to_cartesian(u0)
v = nondivergent_velocity

##
panel_ids = get_panel_ids(panel_model)
degree = 2*(p_fe + 1)

Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,degree)

Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
P = TrialFESpace(Q)

# hard code RT space as order 1 -- for velocity
V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,1); conformity=:HDiv)
U = TrialFESpace(V)

u_cf = panelwise_cellfield(u,Ω_panel,panel_ids)
# _a(u,v) = ∫( u*v )dΩ
# _l(v) = ∫( u_cf*v )dΩ
# op = AffineFEOperator(_a,_l,P,Q)
# uh0 = solve(LUSolver(),op)
uh0 = interpolate(u_cf, P)

meas_cf = CellField(sqrtg,Ω_panel)

# supg stabilisation parameter
_dx = dx(nc(panel_model))
_dt = _dx*CFL/p_fe
dt = floor(_dt,sigdigits=1)

τ = 0.5*dt

function get_velocity(t)
  vecX(XYZ) = v(t)(XYZ)
  vX = panel_to_cartesian(tangent_vec(vecX))
  v_contr_cf =  panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
  # return v_contr_cf
  interpolate(v_contr_cf,U)
end

a_mass_Ω(dtu,v) = ∫( (dtu*v)*meas_cf )dΩ
a_mass_s(t,dtu,v) = ∫( (dtu*(get_velocity(t)⋅∇(v)))*meas_cf )dΩ
a_Ω(t,u,v) = ∫( ((get_velocity(t)⋅∇(u))*v )*meas_cf )dΩ
a_s(t,u,v) =  ∫( ((get_velocity(t)⋅∇(u))*(get_velocity(t)⋅∇(v)) )*meas_cf )dΩ

a_mass(t,dtu,v) = a_mass_Ω(dtu,v) + τ*a_mass_s(t,dtu,v)
res(t,u,v) =  a_Ω(t,u,v) + τ*a_s(t,u,v)
jac(t,u,du,v) = a_Ω(t,du,v) + τ*a_s(t,du,v)
jac_t(t,u,dtu,v) = a_mass_Ω(dtu,v) + τ*a_mass_s(t,dtu,v)
opT = TransientSemilinearFEOperator(a_mass, res, (jac,jac_t), P, Q)

# solve with SSP RK 3
t0, tF = 0.0, 5.0

solver = RungeKutta(LUSolver(), LUSolver(), dt, :EXRK_SSP_3_3)
solT = solve(solver, opT, t0, tF, uh0)

covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)

cell_geo_map = geo_map_func(panel_ids)


writevtk(Ω_panel,dir*"/solT_0.vtu", cellfields=["uh"=>uh0,"v"=>covarient_basis_cf⋅ get_velocity(0.0)],append=false,geo_map=cell_geo_map)


## store errors
ts = Float64[]
Es = Float64[]

push!(ts,0.0)
push!(Es,0.0)

for (t,uh) in solT

  eu = l2((uh-uh0)*meas_cf,dΩ)

  println(t, "; ", eu)

  push!(ts,t)
  push!(Es,eu)

  writevtk(Ω_panel,dir*"/solT_$t.vtu", cellfields=["uh"=>uh,"v"=>covarient_basis_cf⋅ get_velocity(t)],append=false,geo_map=cell_geo_map)

end

_make_pvd_distributed(dir,"solT",1)

output = @strdict ts Es
safesave(datadir(dir, ("advection_errors.jld2")), output)
