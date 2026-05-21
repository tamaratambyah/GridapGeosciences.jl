
using DrWatson
using DataFrames
using Plots
using LaTeXStrings

dir = datadir("grad_grad")
df = collect_results(dir)
lvl = df[(df.p_fe .==1), :n]

ts_am_grad = df[(df.p_fe .==1), :t_ambient]
ts_pm_grad = df[(df.p_fe .==1), :t_panel]

alcs_am_grad = df[(df.p_fe .==1), :alloc_ambient]
alcs_pm_grad = df[(df.p_fe .==1), :alloc_panel]


dir = datadir("coriolis")
df = collect_results(dir)
lvl = df[(df.p_fe .==1), :n]

ts_am_cor = df[(df.p_fe .==1), :t_ambient]
ts_pm_cor = df[(df.p_fe .==1), :t_panel]

alcs_am_cor = df[(df.p_fe .==1), :alloc_ambient]
alcs_pm_cor = df[(df.p_fe .==1), :alloc_panel]


_linestyle = [:solid, :dash, :dot, :dashdotdot, :dashdot, :solid, :dash,]
_markers= [:circle, :rect, :star5, :diamond,  :cross, :hexagon]
_colors = palette(:tab10)

markers= [:circle :rect  :diamond ]
markersize = [6 7 6]
linestyle = [:solid :dash]
colors = [_colors[1] _colors[2] _colors[3]]

xplot = lvl
default(; fontfamily="Computer Modern");
plot(xplot,[ts_am_grad./ts_pm_grad ts_am_cor./ts_pm_cor  ],
    lw=2,
    marker=markers,
    ms = markersize,
    color=colors,
    label=["grad-grad" "coriolis"])
plot!(shape=:auto,
    # xaxis=:log2,
    # yaxis=:log,
    xlabel=latexstring("\$ \\ell \$"),
    ylabel="Time Extrinsic/Instrinsic",
    xtickfontsize=11,ytickfontsize=11,
    xguidefontsize=12,yguidefontsize=12,
    legendfontsize=10,
    legend=:topleft,
    framestyle = :box,
    # guidefontfamily=font(20,"Times Roman"),
    ylimits=(0,6)
    )

xl = map(x->string(Int((x))),xplot)
plot!(xticks = (xplot, xl))
# ys = [1, 1e-3, 1e-6, 1e-9]
# plot!(yticks = (ys))

plot!(show=true)
savefig(plotsdir("benchmark_comparison_time.pdf"))


#########

plot(xplot,[alcs_am_grad./alcs_pm_grad alcs_am_cor./alcs_pm_cor],
    lw=2,
    marker=markers,
    ms = markersize,
    color=colors,
    label=["grad-grad" "coriolis"])
plot!(shape=:auto,
    # xaxis=:log2,
    # yaxis=:log,
    xlabel=latexstring("\$ \\ell \$"),
    ylabel="Allocation Extrinsic/Instrinsic",
    xtickfontsize=11,ytickfontsize=11,
    xguidefontsize=12,yguidefontsize=12,
    legendfontsize=10,
    legend=:topleft,
    framestyle = :box,
    # guidefontfamily=font(20,"Times Roman"),
    # ylimits=(0,22)
    )

xl = map(x->string(Int((x))),xplot)
plot!(xticks = (xplot, xl))
# ys = [1, 1e-3, 1e-6, 1e-9]
# plot!(yticks = (ys))

plot!(show=true)
savefig(plotsdir("benchmark_comparison_allocation.pdf"))
