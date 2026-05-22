using Gridap
using Gridap.Geometry, Gridap.CellData, Gridap.Arrays, Gridap.FESpaces
using Gridap.ReferenceFEs, Gridap.Fields, Gridap.Geometry
using GridapGeosciences
using BenchmarkTools
using GFlops
using DrWatson

include("helper_funcs.jl")
include("ambient_funcs.jl")
include("panel_funcs.jl")


################################################################################
#### Run benchmark
################################################################################
dir = datadir("gradgrad")
degree = 4
orders = collect(1:3)
for order in orders
  benchmark_ambient(order,degree,dir)
  benchmark_reference_panel(order,degree,dir)
end

################################################################################
#### Collect results
################################################################################
using DataFrames
df = collect_results(dir)

ambient = df[df.model.=="ambient",:]
panel = df[df.model.=="panel",:]
# orders = ambient[:,:order]

################################################################################
#### Plot results
################################################################################
using Plots
markers= [:circle :rect  :diamond ]
markersize = [6 7 6]
colors = palette(:tab10)
default(; fontfamily="Computer Modern");

data = [:ops,:t,:flops]
labs = ["Operations", "Time", "Flops"]

plot()
for (i,(sym,lab)) in enumerate(zip(data,labs))
  plot!(orders, ambient[:,sym]./panel[:,sym],
        lw=2,marker=markers[i],ms=markersize[i],color=colors[i],
        label=lab )
end
plot!(show=true)
plot!(shape=:auto,
      xlabel="Order",
      ylabel="Extrinsic/Instrinsic",
      xtickfontsize=11,ytickfontsize=11,
      xguidefontsize=12,yguidefontsize=12,
      legendfontsize=10,
      legend=:bottomright,
      framestyle = :box,
      ylimits=(0,10),
      xticks = (orders, orders)
      )
savefig(plotsdir("benchmark.pdf"))
