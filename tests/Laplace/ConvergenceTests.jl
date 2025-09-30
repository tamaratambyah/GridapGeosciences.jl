using DrWatson
using Gridap
using GridapGeosciences
using Plots, LaTeXStrings

include("Helmholtz.jl")
include("LaplaceBeltrami.jl")
include("analytic_funcs.jl")

n_ref_lvls = 4
return_vtk = false
dir = datadir("LaplaceTests")
(return_vtk && !isdir(dir)) && mkdir(dir)
!isdir(plotsdir()) && mkdir(plotsdir())

### Primial formulation
helmholtz_convergence_test(analytic_funcs,n_ref_lvls,return_vtk)
helmholtz_convergence_test(mapped_funcs,n_ref_lvls,return_vtk)
helmholtz_convergence_test(williamson_funcs,n_ref_lvls,return_vtk)

### Mixed problem
mixed_helmholtz_convergence_test(analytic_funcs,n_ref_lvls,return_vtk)
mixed_helmholtz_convergence_test(mapped_funcs,n_ref_lvls,return_vtk)
mixed_helmholtz_convergence_test(williamson_funcs,n_ref_lvls,return_vtk)

### Laplace-Beltrami - requires zeromean functions
laplace_beltrami_convergence_test(analytic_funcs,n_ref_lvls,return_vtk)
laplace_beltrami_convergence_test(mapped_funcs,n_ref_lvls,return_vtk)
laplace_beltrami_convergence_test(williamson_funcs,n_ref_lvls,return_vtk)
