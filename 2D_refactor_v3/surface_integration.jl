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

include("surface_metric.jl")


# _model = cube_model_3D
_model = ref_ref_ref_model
manifold_grid = CubedSphereGrid(_model)
grid = get_grid(_model) ### underlying coarse grid``
panel_ids = get_panel_ids(_model)

order = 4

manifold_model = UnstructuredDiscreteModel(manifold_grid.parametric_grid)
Ω = Triangulation(manifold_model)
dΩ = Measure(Ω,order)

f = CellField(1,Ω)
sqrt_det_g = CellField(sqrt_det_func,Ω)



quad = CellQuadrature(Ω,order) #dΩ.quad
# trian_f = get_triangulation(f)
# trian_x = get_triangulation(quad)
# trian_sqrt_g = get_triangulation(sqrt_det_g)

# @Gridap.Helpers.check is_change_possible(trian_f,trian_x)
# @Gridap.Helpers.check is_change_possible(trian_sqrt_g,trian_x)


b = change_domain(f,quad.trian,quad.data_domain_style)
sqrt_g = change_domain(sqrt_det_g,quad.trian,quad.data_domain_style)

x = get_cell_points(quad)

bx = b(x)
sqrt_g_x = sqrt_g(x)

# cell_map = get_ambient_cell_map(manifold_grid)
# cell_Jt = lazy_map(∇,cell_map)
# cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)

cell_to_shapefuns = get_cell_shapefuns(manifold_grid.parametric_grid)
cell_to_coords = get_cell_coordinates(manifold_grid)



_cell_map = lazy_map(linear_combination,cell_to_coords,cell_to_shapefuns)

_cell_Jt = lazy_map(∇,_cell_map)
_cell_Jtx = lazy_map(evaluate,_cell_Jt,quad.cell_point)

for i in 1:6
  dd = evaluate(get_cell_map(manifold_grid)[i],quad.cell_point[i]) .== (evaluate(_cell_map[i],quad.cell_point[i]))
  println(sum(dd))
end

## compute the integral by hand
weights = collect1d(quad.cell_weight)
jtx = collect1d(_cell_Jtx)
sgx = collect1d(sqrt_g_x)
_bx = collect1d(bx)
z = 0.0

for j in 1:num_cells(_model)
  aq = _bx[j]
  jq = jtx[j]
  w = weights[j]
  d = sgx[j]
  @inbounds for i in eachindex(aq)
    z+=(aq[i]*w[i]*(Gridap.TensorValues.meas(jq[i]))*d[i]) # multiply by sqrt(det(g))
  end

end
z
4*π*r^2
