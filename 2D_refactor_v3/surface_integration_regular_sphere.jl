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
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools

a = 1.0
r = a*sqrt(3.0)

E(x) = r^2 * 1.0

function G(x)
  u,v = x
  r^2*(cos(u))^2
end

F(x) = 0.0
metric(x) = TensorValue{2,2}(E(x),F(x),F(x),G(x))

order = 4

parametric_model = CartesianDiscreteModel((-π,π,-π/2,π/2),(18,18))
Ω = Triangulation(parametric_model)
dΩ = Measure(Ω,order)

f = CellField(1.0,Ω)
_metric = CellField(metric,Ω)


quad = CellQuadrature(Ω,order)

b = change_domain(f,quad.trian,quad.data_domain_style)
g = change_domain(_metric,quad.trian,quad.data_domain_style)

x = get_cell_points(quad)

bx = b(x)
gx = g(x)

cell_map = get_cell_map(quad.trian)
cell_Jt = lazy_map(∇,cell_map)
cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)

## compute the integral by hand
weights = collect1d(quad.cell_weight)
jtx = collect1d(cell_Jtx)
sgx = map(x-> sqrt.(meas.(x)), gx)
_bx = collect1d(bx)
z = 0.0

for j in 1:num_cells(parametric_model)
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
