using Gridap
using Plots
using LaTeXStrings

include("../src/initialise.jl")


function area_test(manifold_model,order::Int)
  computed_area = get_surface_area(manifold_model,order)
  _true_area = true_area(get_manifold_name(manifold_model))
  area_diff = abs(computed_area - _true_area)

  println("\nTrue area: $(_true_area)\nComputed area: $(computed_area)\nDifference: $(area_diff)")

  @test area_diff < 1e-8
end

################################################################################
#### test surface area of cube
################################################################################
## test the area of cube  -- apply 1 level of refinement
_manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(1),cube)
manifold_model = Adaptivity.refine(_manifold_model)

area_test(manifold_model,2)
area_test(manifold_model,4)

################################################################################
#### surface area of sphere
################################################################################
# ## test the area of cubed sphere -- apply 1 level of refinement
_manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
manifold_model = Adaptivity.refine(_manifold_model)
area_test(manifold_model,2)
area_test(manifold_model,4)
area_test(manifold_model,6)
area_test(manifold_model,8)

################################################################################
#### surface area convergence
################################################################################
nc = [1,2,4,8]  # sqrt num cells per panel
qdegrees = [2,4,6,8] # quadrature orders

# set up plot attributes
markers = [:circle, :rect, :diamond, :utriangle, :cross, :xcross]
_colors = palette(:tab10)

### New cubed sphre model: test convergence over series of refined models
manifold_model = ManifoldDiscreteModel(coarse_cube_model_3D(π/4),cubedsphere)
models = get_refined_models(manifold_model,3)


## compute the convergence rates
dx =   ( sqrt.( 4*π*r^2 ./ (6*nc.^2) ) ) # average width of cell on sphere
gg = [1e-3dx.^2, 5e-6dx.^4, 5e-9dx.^6, 5e-12dx.^8] # some manipulation to get the plot to look nice
leginf = [latexstring("\$ (\\Delta x)^2 \$"), latexstring("\$ (\\Delta x)^4 \$"),
latexstring("\$ (\\Delta x)^6 \$"), latexstring("\$ (\\Delta x)^8 \$")]


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
  lw=4,ms=6,
  c=_colors[i],
  markershape=markers[i],
      label = "q_degree: $degree"
      )

  plot!(0:length(dx)-1, gg[i],
      lw=2,ms=6,
      c=_colors[i],
      label=leginf[i],
      linestyle=:dot,
      yscale=:log10)
end
plot!(yscale=:log10,framestyle=:box,
# title = "surface area of cubed sphere",
xlabel=L"n",
ylabel=L"|a_h - 4πr^2|/4πr^2"
)
# plot!(xlabel="refinement level",
#       ylabel=L"|a_h - 4πr^2|/4πr^2",)
plot!(show=true,legend=:bottomleft,legend_column=2,background_color_legend = nothing)
plot!(xtickfontsize=12,ytickfontsize=12,
legendfontsize=8,guidefontsize=18)
savefig(plotsdir()*"/surface_area_comparison_error")
