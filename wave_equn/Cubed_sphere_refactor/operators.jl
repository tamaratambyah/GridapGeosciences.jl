using DrWatson
using Gridap
using GridapGeosciences
using Test

### First fundamental form
x = Point(0.5,-1,1)
θϕ = xyz2θϕ(x)
g = GridapGeosciences.G_unit_sphere(θϕ)

θ,ϕ = θϕ
@test g[1,1] == (cos(ϕ))^2
@test g[2,1] == g[1,2] == 0.0
@test g[2,2] == 1


test = GridapGeosciences.J_unit_sphere(x)



# gradient
u(x) = x[1] + x[2] + x[3]
