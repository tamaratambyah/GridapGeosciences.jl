using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../src/initialise.jl")

manifold_model = ManifoldDiscreteModel(cube_model_3D,cubedsphere)
Ω = Triangulation(manifold_model)
pts = get_cell_points(Ω)

x = Point(0.5,0.5)

# m = Metric(cubedsphere,Ω)
m = Metric(manifold_model)

## function auto diff
f(x) = x[1]
g(x) = VectorValue(x[1],x[2])

surf_grad = surface_gradient(1.0,m)
surf_grad(x)

surf_grad = surface_gradient(f,m)
surf_grad(x)

surf_div = surface_divergence(g,m)
surf_div(x)

# surf_lap = surface_laplacian(f,m)
# surf_lap(x)


### CellField diff
f_cf = CellField(f,Ω)
g_cf = CellField(g,Ω)

surf_grad = surface_gradient(f_cf,m)
surf_grad(pts)./1

surf_div = surface_divergence(g_cf,m)
surf_div(pts)./1

surf_lap = surface_laplacian(f_cf,m)
surf_lap(pts)./1
