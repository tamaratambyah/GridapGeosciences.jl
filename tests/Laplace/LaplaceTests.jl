using DrWatson
using Gridap
using GridapGeosciences
using Plots, LaTeXStrings

include("Helmholtz.jl")
include("LaplaceBeltrami.jl")
include("analytic_funcs.jl")

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
