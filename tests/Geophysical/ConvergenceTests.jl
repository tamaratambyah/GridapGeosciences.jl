using DrWatson
using Gridap
using GridapGeosciences
using Plots, LaTeXStrings

include("Williamson_functions_v2.jl")

include("wave_equation.jl")
include("shallow_water.jl")
include("nonlinear_shallow_water.jl")

n_ref_lvls = 4

### Williamson 2 tests
williamson2_convergence_test(wave_errors,n_ref_lvls)
williamson2_convergence_test(linear_shallow_water_errors,n_ref_lvls)
williamson2_convergence_test(nonlinear_shallow_water_errors,n_ref_lvls)


### linearised shallow water with arbitary functions
depth(XYZ) = 1.0 + 0.1*exp(-( XYZ[2]^2 + XYZ[3]^2 ) )
velocity(XYZ) = VectorValue(-XYZ[2],XYZ[1],0.0)
coriolis(XYZ) = 2.0

h = panel_to_cartesian(depth)
vecX = velocity
vX = panel_to_cartesian(tangent_vec(vecX))
f = panel_to_cartesian(coriolis)

linear_shallow_water_convergence_test(n_ref_lvls,h,vX,f,true)
