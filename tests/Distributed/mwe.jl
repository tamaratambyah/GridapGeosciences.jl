using DrWatson
using Gridap
using GridapDistributed
using GridapP4est
using GridapGeosciences

using MPI
using PartitionedArrays

include("../convergence_tools.jl")

nprocs = 2
ranks = with_debug() do distribute
  distribute(LinearIndices((nprocs,)))
end

## initial conditions
vecX(XYZ) = zero(XYZ)
depth(XYZ) = 1.0 + 0.01*exp(-5*((1-XYZ[1])^2+(0-XYZ[2])^2+(0-XYZ[3])^2))


h = panel_to_cartesian(depth)
vX = panel_to_cartesian(tangent_vec(vecX))

models  = get_distributed_refined_models(ranks,nprocs,2)
panel_model = models[1]
p_fe = 1

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
h_cf = panelwise_cellfield(h,Ω_panel,panel_ids)
vec_contra_cf = panelwise_cellfield(contra_v(vX),Ω_panel,panel_ids)
xh0 = interpolate([vec_contra_cf,h_cf],X)


## transient weak form
metric_cf = panelwise_cellfield(metric,Ω_panel,panel_ids)
meas_cf = panelwise_cellfield(sqrtg,Ω_panel,panel_ids)
grad_meas_cf = panelwise_cellfield(grad_meas,Ω_panel,panel_ids)

mass(t, (dtu,dtp), (v,q)) = ∫( (v⋅ (metric_cf⋅ dtu) )*meas_cf )dΩ + ∫( (q*dtp)*meas_cf )dΩ
res(t,(u,p),(v,q)) =  ∫( q*(u⋅grad_meas_cf + meas_cf*(∇⋅u) )  )dΩ - ∫( p*(v⋅grad_meas_cf + meas_cf*(∇⋅v) ) )dΩ
jac(t,(u,p),(du,dp),(v,q)) =  res(t,(du,dp),(v,q))
jac_t(t,(u,p),(dtu,dtp),(v,q)) =  mass(t, (dtu,dtp), (v,q))

opT = TransientSemilinearFEOperator(mass,res,(jac,jac_t),X,Y, constant_mass=true)

# transient parameters
t0 = 0.0
dt = 0.01
tF = 1*dt

# solve with SSP RK 3
solver = RungeKutta(LUSolver(), LUSolver(), dt, :EXRK_SSP_3_3)
solT = solve(solver, opT, t0, tF, xh0)

uh = xh0[1]
ph = xh0[2]
for (t, xh) in solT
  uh,ph = xh
  i_am_main(ranks) && println("t = ", t)
end

covarient_basis_cf = panelwise_cellfield(covarient_basis,Ω_panel,panel_ids)
cell_geo_map = geo_map_func(Ω_panel)
dir = datadir("TransientTest")
writevtk(Ω_panel,dir*"/solT_test" * ".vtu", cellfields= ["uh"=>covarient_basis_cf⋅uh, "ph"=>ph],append=false,geo_map=cell_geo_map)


## extract free values
p_vec =  ph.metadata.free_values
u_vec = uh.metadata.free_values

## save free values to file


## load free values from file
# p_free_vals =
# u_free_vals =

## interpoalte onto space
# ph_load = interpolate(p_free_vals, P)
# uh_load = interpolate(u_free_vals,U)
# xh_load = interpolate([uh_load,ph_load],X) ## mutlifield
