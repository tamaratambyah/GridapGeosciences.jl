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
include("darcy_mass_conservation.jl")

n_ref_lvls = 4
ranks = [true]
nprocs = 1
ps = [1,2,3]

### compute surface of sphere
# surface_area_convergence_test(n_ref_lvls)




## scalar tests
for dict in [analytic_funcs] #williamson_funcs ]
  ### interpolation into FE space
  interpolation_convergence_test(ranks,nprocs,dict,n_ref_lvls,ps)

  ### compute sgrad
  sgrad_convergence_test(ranks,nprocs,dict,n_ref_lvls,ps)


end

## vector tests
for dict in [ambient_vecs] #williamson_vec ]
  ## vector projection
  vector_proj_convergence_test(ranks,nprocs,dict,n_ref_lvls,ps)

  ## perp operator
  vector_perp_convergence_test(ranks,nprocs,dict,n_ref_lvls,ps)
end


### mass conversation with scalar fields
mass_conservation_convergence_test(ranks,nprocs,analytic_funcs,true,n_ref_lvls,[1])
mass_conservation_convergence_test(ranks,nprocs,williamson_funcs,true,n_ref_lvls,[1])

### mass conversation with vector fields
mass_conservation_convergence_test(ranks,nprocs,ambient_vecs,false,n_ref_lvls,[1])
mass_conservation_convergence_test(ranks,nprocs,williamson_vec,false,n_ref_lvls,[1])
