using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.Fields
using Gridap.TensorValues
using Gridap.Helpers
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools
include("helpers.jl")
include("maps/metric_maps.jl")
include("surface_metric_and_op/metric_info.jl")
include("surface_metric_and_op/quadrature.jl")
include("surface_metric_and_op/cubedsphere_metric.jl")


r = sqrt(3.0)

order = 4

parametric_model = CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(8,8))
Ω = Triangulation(parametric_model)
dΩ = Measure(Ω,order)

f = CellField(1.0,Ω)
metric = CellField(metric_func,Ω)

quad = CellQuadrature(Ω,order)

b = change_domain(f,quad.trian,quad.data_domain_style)
g = change_domain(metric,quad.trian,quad.data_domain_style)

x = get_cell_points(quad)

bx = b(x)
gx = g(x)

cell_map = get_cell_map(quad.trian)
cell_Jt = lazy_map(∇,cell_map)
cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)

## compute the integral by hand
weights = collect1d(quad.cell_weight)
jtx = collect1d(cell_Jtx)
gx_meas = lazy_map(MetricMeasure(),gx)
_bx = lazy_map(LazyMult(), bx, gx_meas)# multiply by sqrt(det(g))
z = 0.0

for j in 1:num_cells(parametric_model)
  aq = _bx[j]
  jq = jtx[j]
  w = weights[j]
  d = gx_meas[j]
  @inbounds for i in eachindex(aq)
    z+=(aq[i]*w[i]*(Gridap.TensorValues.meas(jq[i]))*d[i] )
  end

end
z
6*z # assume each face of cube has same area
4*π*r^2

s_quad = SurfaceQuadrature(metric_func,quad)
6*sum(integrate(1.0,s_quad))
4*π*r^2

_s_quad = SurfaceQuadrature(metric,quad)
6*sum(integrate(1.0,_s_quad))


m = MetricInfo(metric_func,Ω)
_s_quad = SurfaceQuadrature(m,quad)
6*sum(integrate(1.0,_s_quad))

DΩg = Measure(s_quad)
6*sum( integrate(1.0,DΩg))

dΩg = Measure(m,Ω,order)
6*sum( integrate(1.0,dΩg))

6*sum( ∫(1  )dΩg )
