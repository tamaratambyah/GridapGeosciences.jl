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

_dirs = ["gadi/LaplaceBeltramiConvergence/func_sin",
        "gadi/LaplaceBeltramiConvergence/func_XYZ",
        "gadi/HelmholtzConvergence/func_sin",
        "gadi/HelmholtzConvergence/func_XYZ",
        "gadi/WaveConvergence/func_z1",
        "gadi/AdvectionSUPGConvergence",
        # "gadi/TransientAdvectionSUPGConvergence_Octree",
        "gadi/ShallowWaterConvergence/func_z1",
         "gadi/TransientShallowWaterConvergence"
]
# dir = datadir()
# dir = datadir("gadi/AdvectionDGConvergence")
# dir = datadir()
# dir = datadir("gadi/TransientAdvectionDGUpwinding")
for d in _dirs
  dir = datadir(d)
  files = filter(x->endswith(x, ".jld2"), readdir((dir)))
  plot_convergence_from_saved(dir,"convergence",["u","p","n"])
end


# ### analysis: Laplace
using DrWatson
using DataFrames
include("convergence_tools.jl")
dir = datadir("LaplaceBeltramiConvergence_3D/func_MXYZ/convergence")
# dir = datadir("LaplaceBeltramiConvergence_3D/func_Msin/convergence")
df = collect_results(dir)

ps = unique(df.p_fe)
plot()
for p in ps
  errors = df[(df.p_fe .== p ),:e]
  dxs = df[(df.p_fe .== p ),:dxx]
  ns = df[(df.p_fe .== p ),:n]

  slope = convergence_rate(dxs,errors)
  plot_convergence(errors,ns,dxs,slope;leginf=["u"],colors=[palette(:tab10)[p],palette(:tab10)[p] ] )
end
plot!(show=true)
savefig(dir*"/convergence_laplace_3D")



# ### analysis: Wave equation
using DrWatson
using DataFrames
# using GridapGeosciences
include("convergence_tools.jl")
dir = datadir("WaveConvergence_3D/func_z1/convergence")
# dir = datadir("LinearisedShallowWater_3D/func_z1/convergence")
df = collect_results(dir)

ps = unique(df.p_fe)

plot()
for p in ps
  e_u = df[(df.p_fe .== p ),:e_u]
  e_p = df[(df.p_fe .== p ),:e_p]

  dxs = df[(df.p_fe .== p ),:dxx]
  ns = df[(df.p_fe .== p ),:n]

  slope_u = convergence_rate(dxs,e_u)
  slope_p = convergence_rate(dxs,e_p)
  errors = [e_u;e_p]
  plot_convergence(errors,ns,dxs,slope_u;leginf=["u","p"],colors=[palette(:tab10)[p],palette(:tab10)[p] ] )
end

plot!(show=true)
savefig(dir*"/wave_equation_3D")



using DrWatson
using DataFrames
include("convergence_tools.jl")
dir = datadir("LinearBoussineseqConvergence/convergence")
df = collect_results(dir)

ps = unique(df.p_fe)

plot()
for p in ps
  e_u = df[(df.p_fe .== p ),:e_u]
  e_p = df[(df.p_fe .== p ),:e_p]
  e_b = df[(df.p_fe .== p ),:e_b]

  dxs = df[(df.p_fe .== p ),:dxx]
  ns = df[(df.p_fe .== p ),:n]

  slope_u = convergence_rate(dxs,e_u)
  slope_p = convergence_rate(dxs,e_p)
  errors = [e_u;e_p;e_b]
  plot_convergence(errors,ns,dxs,slope_u;leginf=["u","p","b"],
    colors=[palette(:tab10)[p],palette(:tab10)[p],palette(:tab10)[p] ] )
end
plot!(show=true)
savefig(dir*"/linear_boussineseq_3D")


dir = datadir("gadi/sol_p1_nref_h3_nref_v3")
# τ = (1/7.292e-5)/(3600*24)
_make_pvd_distributed(dir,"solT",1.0)

include("plot_tools.jl")
dir = datadir("gadi/TransientShallowWater_W5_octree_supg/sol_p1_nref6")
n = 100
plotName = "SW"
cf = "ph"
plot_latlon(dir,n,plotName,cf)
