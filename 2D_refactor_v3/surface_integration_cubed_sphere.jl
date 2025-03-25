using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.Fields
using Gridap.Helpers
using Gridap.TensorValues
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools

include("initialise.jl")
include("maps/metric_maps.jl")
include("surface_metric.jl")


# _model = cube_model_3D
_model = ref_ref_model
manifold_grid = CubedSphereGrid(_model)

order = 6

manifold_model = UnstructuredDiscreteModel(manifold_grid.parametric_grid)
Ω = Triangulation(manifold_model)
dΩ = Measure(Ω,order)

f = CellField(1.0,Ω)
g = CellField(metric_func,Ω)

quad = CellQuadrature(Ω,order)

b = change_domain(f,quad.trian,quad.data_domain_style)

x = get_cell_points(quad)
bx = b(x)

_g = change_domain(g,quad.trian,quad.data_domain_style)
x = get_cell_points(quad)
gx = _g(x)
gx_meas = lazy_map(MetricMeasure(),gx) # get sqrt(meas(g)) for the area element



node_to_coords = get_node_coordinates(manifold_grid.ambient_grid)
cell_to_nodes = get_cell_node_ids(manifold_grid)
lazy_map(Broadcasting(Reindex(node_to_coords)),cell_to_nodes)


cell_to_coords = collect1d( lazy_map(evaluate,get_cell_map(manifold_grid),get_cell_ref_coordinates(manifold_grid)) )

# cell_to_coords = manifold_grid.parametric_cell_coords


type_to_reffes = get_reffes(manifold_grid)
cell_to_type = get_cell_type(manifold_grid)
type_to_shapefuns = map(get_shapefuns, type_to_reffes)
cell_to_shapefuns = expand_cell_data(type_to_shapefuns,cell_to_type)



_cell_map = lazy_map(linear_combination,cell_to_coords,cell_to_shapefuns)



_cell_Jt = lazy_map(∇,_cell_map)
_cell_Jtx = lazy_map(evaluate,_cell_Jt,quad.cell_point)





## compute the integral by hand
weights = collect1d(quad.cell_weight)
jtx = collect1d(_cell_Jtx)
bgx = lazy_map(LazyMult(), bx,  gx_meas) # multiply by sqrt(det(g))
z = 0.0

for j in 1:num_cells(manifold_model)
  aq = bgx[j]
  jq = jtx[j]
  w = weights[j]
  @inbounds for i in eachindex(aq)
    z+=(aq[i]*w[i]*(Gridap.TensorValues.meas(jq[i])) )
  end

end
z

4*π*r^2
