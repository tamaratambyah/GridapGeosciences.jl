
using DrWatson
using DataFrames
using Plots
using LaTeXStrings

dir = datadir("grad_grad_single_panel")


df = collect_results(dir)
lvl = df[(df.p_fe .==1), :n]

ts_am = df[(df.p_fe .==1), :t_ambient]
ts_pm = df[(df.p_fe .==1), :t_panel]
ts_ref = df[(df.p_fe .==1), :t_panel_ref]

alcs_am = df[(df.p_fe .==1), :alloc_ambient]
alcs_pm = df[(df.p_fe .==1), :alloc_panel]
alcs_ref = df[(df.p_fe .==1), :alloc_panel_ref]


_linestyle = [:solid, :dash, :dot, :dashdotdot, :dashdot, :solid, :dash,]
_markers= [:circle, :rect, :star5, :diamond,  :cross, :hexagon]
_colors = palette(:tab10)

markers= [:circle :rect  :diamond ]
markersize = [6 7 6]
linestyle = [:solid :dash]
colors = [_colors[1] _colors[2] _colors[3]]

xplot = lvl
default(; fontfamily="Computer Modern");
plot(xplot,[ts_am  ts_pm ts_ref],
    lw=2,
    marker=markers,
    ms = markersize,
    color=colors,
    label=["Extrinsic " "Instrinsic" "Reference"])
# plot(xplot,[ts_am./ts_pm ts_pm./ts_ref],
#     lw=2,
#     marker=markers,
#     ms = markersize,
#     color=colors,
#     label=["Extrinsic/Instrinsic" "Instrinsic/Reference"])
plot!(shape=:auto,
    # xaxis=:log2,
    yaxis=:log,
    xlabel=latexstring("\$ \\ell \$"),
    ylabel="Time",
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
savefig(plotsdir("benchmark_"*basename(dir)*"_time.pdf"))


#########

plot(xplot,[alcs_am alcs_pm alcs_ref],
    lw=2,
    marker=markers,
    ms = markersize,
    color=colors,
    label=["Extrinsic" "Instrinsic" "Reference"])
plot!(shape=:auto,
    # xaxis=:log2,
    yaxis=:log,
    xlabel=latexstring("\$ \\ell \$"),
    ylabel="Allocation",
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
savefig(plotsdir("benchmark_"*basename(dir)*"_allocation.pdf"))
