using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.Fields
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools

include("initialise.jl")

_model = cube_model_3D
manifold_grid = CubeGrid(_model)
panel_ids = get_panel_ids(manifold_grid)

Ω = GenericTriangulation(manifold_grid)
dΩ = Measure(Ω,2)

f = CellField(1,Ω)

quad = dΩ.quad
trian_f = get_triangulation(f)
trian_x = get_triangulation(quad)

@Gridap.Helpers.check is_change_possible(trian_f,trian_x)

b = change_domain(f,quad.trian,quad.data_domain_style)
x = get_cell_points(quad)
bx = b(x)


cell_map = get_cell_map(quad.trian)
cell_Jt = lazy_map(∇,cell_map)
cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)

cache = array_cache(cell_Jt)
for i in eachindex(cell_Jt)
  grad = getindex!(cache, cell_Jt, i)
  evaluate(grad,quad.cell_point[i])
end

evaluate(cell_Jt[1],quad.cell_point[1])


lazy_map(IntegrationMap(),bx,quad.cell_weight,cell_Jtx)


# integrate(f,dΩ)
