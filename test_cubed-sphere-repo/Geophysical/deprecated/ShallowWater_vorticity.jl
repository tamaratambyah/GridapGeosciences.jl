"""
solve the non-linear shallow water equations in steady form using manufactured solutions
u + q F^⟂ + ∇ᵧ(Φ) = f₁
φ + ∇ᵧ⋅F = f₁
F = φu
Φ = 0.5(u⋅u) + gᵣφ
q = 1/φ( ∇ᵧ^⟂⋅u  + f )
"""

using DrWatson
using Gridap
using GridapDistributed
using GridapSolvers
using PartitionedArrays
using Gridap.Geometry, Gridap.Adaptivity, Gridap.Helpers, Gridap.Algebra

using GridapGeosciences
using Test

include("../convergence_tools.jl")
include("Williamson2Test.jl")


n_ref_lvls = 3
p_fe = 1
ζ = 0.0
ls = LUSolver()
models  = get_refined_models(n_ref_lvls)

h = panel_to_cartesian(h₀(ζ))
vX = panel_to_cartesian(tangent_vec(u₀(ζ)))
f = panel_to_cartesian(f₀(ζ))
η = panel_to_cartesian(η₀(ζ))


panel_model = models[1]


panel_ids = get_panel_ids(panel_model)
Ω_panel = Triangulation(panel_model)
dΩ = Measure(Ω_panel,2*(p_fe+1))
Ω_error = Triangulation(panel_model)
dΩ_error = Measure(Ω_error,6*p_fe+1)

R = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe+1); conformity=:H1)
H = TrialFESpace(R)

Q = TestFESpace(panel_model, ReferenceFE(lagrangian,Float64,p_fe); conformity=:L2)
P = TrialFESpace(Q)

V = TestFESpace(panel_model, ReferenceFE(raviart_thomas,Float64,p_fe); conformity=:HDiv)
U = TrialFESpace(V)


## initial conditions
covariant_basis_cf = ParametricCellField(covariant_basis,Ω_panel,panel_ids)
u_contra_cf = ParametricCellField(contra_v(vX),Ω_panel,panel_ids)
u_contra_h = interpolate(u_contra_cf,U)
# u_proj_h = covariant_basis_cf ⋅ u_contra_h
u_proj_h = covariant_basis_cf ⋅ u_contra_cf

h_cf = ParametricCellField(h,Ω_panel,panel_ids)
h_h = interpolate(h_cf,P)

cor_cf = ParametricCellField(f,Ω_panel,panel_ids)
gravity = _g

# absolute vorticity
η_cf = ParametricCellField(η,Ω_panel,panel_ids)
η_h = interpolate(η_cf,H)

# mectrics required in weak forms
detg_cf = ParametricCellField(detg,Ω_panel,panel_ids)
metric_cf = ParametricCellField(metric,Ω_panel,panel_ids)
inv_metric_cf = ParametricCellField(inv_metric,Ω_panel,panel_ids)
meas_cf = ParametricCellField(sqrtg,Ω_panel,panel_ids)
grad_meas_cf = ParametricCellField(grad_meas,Ω_panel,panel_ids)

# vorticity
perp_matrix_cf = ParametricCellField(perp_matrix,Ω_panel,panel_ids)
biformq(q,r) = ∫( q*r*meas_cf  )dΩ
liformq(r) = ∫( cor_cf*r*meas_cf  )dΩ + ∫( (perp_matrix_cf⋅u_contra_cf)⋅∇(r)  )dΩ
op = AffineFEOperator(biformq,liformq,H,R)
qh = solve(ls,op)

e_η = l2((η_cf - qh ),meas_cf,dΩ_error)

dir = datadir("voriticity")
!isdir(dir) && mkdir(dir)
cell_geo_map = latlon_geo_map_func(Ω_panel)
panel_cfs = [η_cf, qh, η_cf - qh]
labels = ["eta","etah","ee"]
cellfields = map((x,y) -> x=>y, labels,panel_cfs)
writevtk(Ω_panel,dir*"/vorticity",cellfields=cellfields,append=false,geo_map=cell_geo_map)
