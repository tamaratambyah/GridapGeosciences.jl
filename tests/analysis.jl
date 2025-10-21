using DrWatson
using JLD2
using Plots
using GridapGeosciences

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


#### shallow water
dir = datadir("gadi/Transient_shallow_water_W5_supg/sol_nref6")
_make_pvd_distributed(dir,"solT",1)
df = load(datadir(dir, ("shallow_water_errors.jld2")))

Es_u = df["Es_u"]
Es_p = df["Es_p"]
plot()
plot!(1:length(Es_u),Es_u,lw=3,label="Eu")
plot!(xlabel="t",ylabel="Eu")

plot()
plot!(1:length(Es_p),Es_p,lw=3,label="Ep")
plot!(xlabel="t",ylabel="Ep")


## gadi convegence results
include("convergence_tools.jl")
# dir = datadir("gadi/LaplaceBeltramiConvergence/func_sin")
# dir = datadir("gadi/HelmholtzConvergence/func_sin")
# dir = datadir("gadi/WaveConvergence/func_z1")
# dir = datadir("gadi/AdvectionSUPGConvergence")
# dir = datadir("gadi/AdvectionDGConvergence")
dir = datadir("gadi/TransientAdvectionSUPGConvergence_Octree")
# dir = datadir("gadi/ShallowWaterConvergence/func_z1")
# dir = datadir("gadi/TransientAdvectionDGUpwinding")
dir = datadir("gadi/TransientShallowWaterConvergence")
files = filter(x->endswith(x, ".jld2"), readdir((dir)))
plot_convergence_from_saved(dir,"convergence",["u","p","n"])
# plot_convergence_from_saved(dir,"convergence",["p","u","n"])




include("plot_tools.jl")
dir = datadir("TransientAdvectionSUPG/sol_p1_nref4")
n = 100
plotName = "Advection"
cf = "uh"
plot_latlon(dir,n,plotName,cf)
