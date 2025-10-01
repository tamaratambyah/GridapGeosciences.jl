using DrWatson
using Gridap
using GridapGeosciences
using Plots, LaTeXStrings

include("Helmholtz.jl")
include("LaplaceBeltrami.jl")
include("analytic_funcs.jl")

n_ref_lvls = 4
ps = [1]#,2,3]
ls = LUSolver()
return_vtk = false
dir = datadir("LaplaceTests")
!isdir(dir) && mkdir(dir)


for dict in [analytic_funcs]#, mapped_funcs, williamson_funcs]
  ### Helmholtz: primial formulation
  helmholtz_convergence_test(dir,dict,n_ref_lvls,ps,ls,return_vtk)

  ### Helmholtz: mixed problem
  mixed_helmholtz_convergence_test(dir,dict,n_ref_lvls,ps,ls,return_vtk)

  ### Laplace-Beltrami - requires zeromean functions
  laplace_beltrami_convergence_test(dir,dict,n_ref_lvls,ps,ls,return_vtk)
end
