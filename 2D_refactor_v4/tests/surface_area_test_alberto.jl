using DrWatson
using GridapGeosciences
using Gridap
using Plots
using LaTeXStrings

include("../src/initialise.jl")


## for Alberto's model
function get_surface_area(model::GridapGeosciences.AnalyticalMapCubedSphereDiscreteModel,degree::Int)
  Ω = Triangulation(model)
  dΩ = Measure(Ω,degree)
  computed_area = sum( integrate(1.0,dΩ))
  return computed_area
end

nc = [1,2,4,8]  # sqrt num cells per panel
qdegrees = [2,4,6,8] # quadrature orders



# set up plot attributes
markers = [:circle, :rect, :diamond, :utriangle, :cross, :xcross]
_colors = palette(:tab10)

### New cubed sphre model: test convergence over series of refined models
manifold_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
models = get_refined_models(manifold_model,3)

plot()
# compute and plot errors in surface area
for i in 1:length(qdegrees)
  degree = qdegrees[i]
  errs = []
  for i in 1:length(models)
    manifold_model = models[i]
    err_h = abs(get_surface_area(manifold_model,degree) - true_area(cubedsphere))/ true_area(cubedsphere)
    push!(errs, err_h )
  end
  plot!(0:length(errs)-1,errs,
  lw=3,ms=6,
  c=_colors[i],
  markershape=markers[i],
      label = "q_degree: $degree"
      )
end



#### Alberto's model
for i in 1:length(qdegrees)
  degree = qdegrees[i]
  errs = []
  for i in 1:length(nc)
    model =  CubedSphereDiscreteModel(nc[i];radius=r)
    err_h = abs(get_surface_area(model,degree) -  true_area(cubedsphere))/ true_area(cubedsphere)
    push!(errs, err_h )
  end
  plot!(0:length(errs)-1,errs,
  lw=3,ms=6,
  c=_colors[i],
  label="",
  markershape=markers[i],
  linestyle=:dash)
end
plot!(yscale=:log10,framestyle=:box,
title = "surface area of cubed sphere",)
plot!(xlabel="refinement level",
      ylabel=L"|a_h - 4πr^2|/4πr^2",)


## compute the convergence rates
dx =   ( sqrt.( 4*π*r^2 ./ (6*nc.^2) ) ) # average width of cell on sphere
gg = [1e-3dx.^2, 5e-6dx.^4, 5e-9dx.^6, 5e-12dx.^8] # some manipulation to get the plot to look nice
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
