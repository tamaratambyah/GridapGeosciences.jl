using DrWatson
using Gridap
using GridapGeosciences
using Plots, LaTeXStrings

include("Helmholtz.jl")
include("LaplaceBeltrami.jl")
include("analytic_funcs.jl")

# analytic functions defined on the chart
analytic_funcs = Dict{Symbol,Any}()
analytic_funcs[:sin] = f_sin
analytic_funcs[:XYZ] = f_XYZ

# mapped functions from cartesian or latlon
mapped_funcs = Dict{Symbol,Any}()
mapped_funcs[:sin] = panel_to_latlon(fθϕ)
mapped_funcs[:XYZ] = panel_to_cartesian(fX)


#### Williamson 2 streamfunction
williamson_funcs = Dict{Symbol,Any}()
williamson_funcs[:z1] = panel_to_latlon(fWilliamson(0))
williamson_funcs[:z2] = panel_to_latlon(fWilliamson(0.05))
williamson_funcs[:z3] = panel_to_latlon(fWilliamson(π/2-0.05))
williamson_funcs[:z4] = panel_to_latlon(fWilliamson(π/2))

n_ref_lvls = 4

### Primial formulation
helmholtz_convergence_test(analytic_funcs,n_ref_lvls,false)
helmholtz_convergence_test(mapped_funcs,n_ref_lvls)
helmholtz_convergence_test(williamson_funcs,n_ref_lvls,true)

### Mixed problem
mixed_helmholtz_convergence_test(analytic_funcs,n_ref_lvls,false)
mixed_helmholtz_convergence_test(mapped_funcs,n_ref_lvls)
mixed_helmholtz_convergence_test(williamson_funcs,n_ref_lvls,true)

### Laplace-Beltrami - requires zeromean functions
laplace_beltrami_convergence_test(analytic_funcs,n_ref_lvls,true)
laplace_beltrami_convergence_test(mapped_funcs,n_ref_lvls,true)
laplace_beltrami_convergence_test(williamson_funcs,n_ref_lvls,true)
