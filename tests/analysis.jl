using DrWatson
using JLD2
using Plots

#### advection
dir = datadir("gadi/Transient_advection_nref6")
df = load(datadir(dir, ("advection_errors.jld2")))
ts = df["ts"]
Es = df["Es"]

plot()
plot!(ts,Es,lw=3)
plot!(xlabel="t",ylabel="L2(u - uh)")
savefig(plotsdir()*"/advection_errors")


### wave equation
dir = datadir("Transient_wave_equation")
df = load(datadir(dir, ("wave_errors.jld2")))
ts = df["ts"]
Es = df["Es"]
ms = df["ms"]
s_divus = df["s_divus"]
divus = df["divus"]


# plot relative error in mass and energy
ms_rel = abs.(ms.-ms[1])./ms[1]
Es_rel = abs.(Es.-Es[1])./Es[1]

plot()
plot!(ts[2:end],ms_rel[2:end],lw=3,label="mass")
plot!(ts[2:end],Es_rel[2:end],lw=3,label="energy")
plot!(yaxis=:log,xlabel="t",ylabel=L"|x_t-x_0|/x_0")
savefig(plotsdir()*"/wave_transient_conservation")

# plot surface divergence
plot()
plot!(ts[2:end],abs.(s_divus[2:end]),lw=3)
plot!(yaxis=:log,xlabel="t",ylabel="|s_div (u)|")
savefig(plotsdir()*"/wave_transient_s_div")

# plot panel divergence
plot()
plot!(ts[2:end],abs.(divus[2:end]),lw=3)
plot!(xlabel="t",ylabel="|div (u)|")
savefig(plotsdir()*"/wave_transient_panel_div")
