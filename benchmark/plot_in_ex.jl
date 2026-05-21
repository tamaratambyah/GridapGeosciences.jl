
using DrWatson
using DataFrames
using Plots
using LaTeXStrings



_linestyle = [:solid, :dash, :dot, :dashdotdot, :dashdot, :solid, :dash,]
_markers= [:circle, :rect, :star5, :diamond,  :cross, :hexagon]
_colors = palette(:tab10)

markers= [:circle :rect  :diamond ]
markersize = [6 7 6]
linestyle = [:solid :dash]
colors = [_colors[1] _colors[2] _colors[3]]


dir = datadir("coriolis/flops")


df = collect_results(dir)
lvl = df[:, :nref]

ambient = df[:, :ambient]
panel = df[:, :panel]
ref = df[:, :ref]



xplot = lvl
default(; fontfamily="Computer Modern");
plot(xplot,[ambient panel ref],
    lw=2,
    marker=markers,
    ms = markersize,
    color=colors,
    label=["Extrinsic " "Instrinsic" "Reference"])

plot!(shape=:auto,
    # xaxis=:log2,
    # yaxis=:log,
    xlabel=latexstring("\$ \\ell \$"),
    ylabel="G "*basename(dir),
    xtickfontsize=11,ytickfontsize=11,
    xguidefontsize=12,yguidefontsize=12,
    legendfontsize=10,
    legend=:topleft,
    framestyle = :box,
    # guidefontfamily=font(20,"Times Roman"),
    # ylimits=(1e-9,0)
    )

xl = map(x->string(Int((x))),xplot)
plot!(xticks = (xplot, xl))
# ys = [1, 1e-3, 1e-6, 1e-9]
# plot!(yticks = (ys))

plot!(show=true)
savefig(plotsdir("benchmark_"*basename(dir)*".pdf"))
