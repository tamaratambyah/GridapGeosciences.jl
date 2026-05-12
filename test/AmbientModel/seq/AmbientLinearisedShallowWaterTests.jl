include("../AmbientLinearisedShallowWater.jl")

## Serial model: 2D
n_ref_lvls = 4
radius = 1.0
models = get_ambient_refined_models(n_ref_lvls,radius)
AmbientLinearisedShallowWaterTests.main(models)



# ### I do not like having this here, but need to think of a better way
# to compare one error result to the intrinsic approach
using Gridap
include("../../Geophysical/Williamson_functions.jl")
a_e = 6.37e6 # m
g = 9.8 # m/2
ω = 7.29e-5 #s^-1
T = 12*24*3600 #s
H_0 = 2.94e4/g #m
u_0 = 2*π*a_e/T #m/s

L = a_e
_τ = 1/ω

_a = a_e/L
_g = g*_τ^2/L
_ω = ω*_τ
_H_0 = H_0/L
_T = T/_τ
_u0 = u_0/L*_τ

dir = @__DIR__
p_fe = 1


ambient_model = models[3]
h = h₀(0.0)
vX = tangent_vec(u₀(0.0))
f = f₀(0.0)
e_u_ambient, e_p_ambient = AmbientLinearisedShallowWaterTests.linear_shallow_water_solver(
  ambient_model,p_fe,dir,
  h,vX,f)


include("../../Geophysical/LinearisedShallowWater.jl")


h = panel_to_cartesian(h₀(0.0))
vX = panel_to_cartesian(tangent_vec(u₀(0.0)))
f = panel_to_cartesian(f₀(0.0))
panel_model = get_parametric_model(ambient_model)
e_u_panel, e_p_panel = LinearisedShallowWaterTests.linear_shallow_water_solver(
  panel_model,p_fe,dir,
  h,vX,f)



e_u_comparison = e_u_ambient - e_u_panel
println(e_u_comparison)
@test e_u_comparison < 1e-12

e_p_comparison = e_p_ambient - e_p_panel
println(e_p_comparison)
@test e_p_comparison < 1e-12
