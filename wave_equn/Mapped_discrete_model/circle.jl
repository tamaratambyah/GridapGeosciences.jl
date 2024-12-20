using DrWatson
using Plots


t = LinRange(0,π,10)

xs = []
ys = []
for i in 1:length(t)
  push!(xs,cos(t[i]))
  push!(ys,sin(t[i]))
end

plot()
plot!(xs,ys)
