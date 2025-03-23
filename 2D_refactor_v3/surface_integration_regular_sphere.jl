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

a = 1.0
r = a*sqrt(3.0)


#### need the sqrt(det(metric)) in the integral
function sqrt_det_func(x)
 u,v = x
 sqrt( r^4*(cos(u))^2 )
end

order = 10

parametric_model = CartesianDiscreteModel((-π,π,-π/2,π/2),(10,10))
Ω = Triangulation(parametric_model)
dΩ = Measure(Ω,order)

f = CellField(1.0,Ω)
sqrt_det_g = CellField(sqrt_det_func,Ω)

quad = CellQuadrature(Ω,order)

b = change_domain(f,quad.trian,quad.data_domain_style)
sqrt_g = change_domain(sqrt_det_g,quad.trian,quad.data_domain_style)

x = get_cell_points(quad)

bx = b(x)
sqrt_g_x = sqrt_g(x)

cell_map = get_cell_map(quad.trian)
cell_Jt = lazy_map(∇,cell_map)
cell_Jtx = lazy_map(evaluate,cell_Jt,quad.cell_point)

## compute the integral by hand
weights = collect1d(quad.cell_weight)
jtx = collect1d(cell_Jtx)
sgx = collect1d(sqrt_g_x)
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
