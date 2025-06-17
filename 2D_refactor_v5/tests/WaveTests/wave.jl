
using Gridap, Gridap.Geometry
using Plots, LaTeXStrings
using LinearAlgebra
using DrWatson
dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


include("wave_helpers.jl")
include("../PoissonTests/poisson_helpers.jl")



p = 2
degree = 2*(p+1)
ns = [2^i for i = 2:6]
dx = 1 ./ ns

radii = [0.1, 0.5, 1.0, 2.0, 4.0]

################################################################################
#### 1D tests: circle
################################################################################
domain = (0.0, 2*π)
partition = (n, )# n)
uex(x) = VectorValue(x[1]*(2π-x[1]))
pex(x) = 1.0 + 0.01*x[1]*(2π-x[1])

errs_u = []
errs_ug = []
errs_p = []
errs_pg = []

for r in collect(radii)
  metric_func(x) = TensorValue{1}(r^2)
  for n in collect(ns)
    println(n)

    eu,ep,eu_g,ep_g = solve_wave_manifold((-π/2,π/2 ,0,2π), (n, ),p,degree,metric_func,uex,pex)

    push!(errs_u,eu)
    push!(errs_ug,eu_g)
    push!(errs_p,ep)
    push!(errs_pg,ep_g)
  end
end

################################################################################
#### 2D tests: sphere
################################################################################
domain = (-π/2,π/2,0,2π)
partition = (n,n)
metric_func(x) = TensorValue{2,2}(radius^2 ,0.0, 0.0, radius^2*(cos(x[1]))^2 )
uex(x) = VectorValue(0.0,x[2]*(2π-x[2]))
pex(x) = 1.0 +  0.01*(π/2+x[1])*(π/2-x[1])
