using DrWatson
using Gridap
using GridapGeosciences
using Plots, LaTeXStrings

include("interpolation.jl")
include("sgradient.jl")
include("surface_area.jl")
include("vector_projection.jl")
include("vector_perp.jl")
include("analytic_funcs.jl")


n_ref_lvls = 4

### compute surface of sphere
surface_area_convergence_test(n_ref_lvls)

### interpolation into FE space
interpolation_convergence_test(analytic_funcs,n_ref_lvls)
interpolation_convergence_test(williamson_funcs,n_ref_lvls)

### Computation of sgrad
sgrad_convergence_test(analytic_funcs,n_ref_lvls)
sgrad_convergence_test(williamson_funcs,n_ref_lvls)

### vector projection
vector_proj_convergence_test(ambient_vecs,n_ref_lvls,false)
vector_proj_convergence_test(williamson_vec,n_ref_lvls,true)

### perp operator
vector_perp_convergence_test(ambient_vecs,n_ref_lvls,false)
vector_perp_convergence_test(williamson_vec,n_ref_lvls,true)
