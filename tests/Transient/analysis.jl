using DrWatson
using JLD2
using Plots

dir = datadir("gadi/Transient_advection_nref6")
df = load(datadir(dir, ("advection_errors.jld2")))
ts = df["ts"]
Es = df["Es"]

plot()
plot!(ts,Es,lw=3)
plot!(xlabel="t",ylabel="L2(u - uh)")
