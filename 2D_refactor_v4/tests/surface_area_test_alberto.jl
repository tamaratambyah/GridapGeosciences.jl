using DrWatson
using GridapGeosciences
using Gridap
using Plots
using LaTeXStrings

include("../src/initialise.jl")


## for Alberto's model
function get_surface_area(model::GridapGeosciences.AnalyticalMapCubedSphereDiscreteModel,order::Int)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,order)
  computed_area = sum( integrate(1.0,dΩ))
  return computed_area
end

nc = [1,2,4,8]  # sqrt num cells per panel
orders = [1,2,3,4,5,6] # quadrature orders



# set up plot attributes
markers = [:circle, :rect, :diamond, :utriangle, :cross, :xcross]
_colors = palette(:tab10)

### New cubed sphre model: test convergence over series of refined models
manifold_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
models = get_refined_models(manifold_model,3)

plot()
# compute and plot errors in surface area
for order in orders
  errs = []
  for i in 1:length(models)
    manifold_model = models[i]
    err_h = abs(get_surface_area(manifold_model,order) - true_area(cubedsphere))
    push!(errs, err_h )
  end
  plot!(0:length(errs)-1,errs,
  lw=3,ms=6,
  c=_colors[order],
  markershape=markers[order],
      label = "q_order: $order"
      )
end



#### Alberto's model
for order in orders
  errs = []
  for i in 1:length(nc)
    model =  CubedSphereDiscreteModel(nc[i];radius=r)
    err_h = abs(get_surface_area(model,order) -  true_area(cubedsphere))
    push!(errs, err_h )
  end
  plot!(0:length(errs)-1,errs,
  lw=3,ms=6,
  c=_colors[order],
  label="",
  markershape=markers[order],
  linestyle=:dash)
end
plot!(yscale=:log10,framestyle=:box,
title = "surface area of cubed sphere",)
plot!(xlabel="refinement level",
      ylabel=L"|a_h - 4πr^2|",)


## compute the convergence rates
dx =   ( sqrt.( 4*π*r^2 ./ (6*nc.^2) ) ) # average width of cell on sphere
gg = [dx.^1, 1e-2dx.^2, 1e-2dx.^3, 1e-2dx.^4, 5e-4dx.^4, 5e-6dx.^6] # some manipulation to get the plot to look nice
for i in 1:length(gg)
  plot!(0:length(dx)-1, gg[i],
  lw=2,ms=6,
  c=_colors[i],
  label="",
  linestyle=:dot,
  yscale=:log10)
end
plot!(show=true)

savefig(plotsdir()*"/surface_area_comparison_error")
