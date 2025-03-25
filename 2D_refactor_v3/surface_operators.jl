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
using Gridap.CellData
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools
include("helpers.jl")
include("surface_metric.jl")

parametric_model = UnstructuredDiscreteModel(CartesianDiscreteModel((-π,π,-π/2,π/2),(8,8)))
Ω = Triangulation(parametric_model)
pts = get_cell_points(Ω)



struct MetricInfo{A<:CellField,B<:CellField,C<:CellField} <: CellDatum
  metric::A
  sq_meas::B
  inv_metric::C
end

function MetricInfo(metric::CellField,sq_meas::CellField,inv_metric::CellField)
  A = typeof(metric)
  B = typeof(sq_meas)
  C = typeof(inv_metric)
  MetricInfo{A,B,C}(metric,sq_meas,inv_metric)
end

function MetricInfo(metric_func::Function,Ω::Triangulation)
  sq_meas_func(x) = sqrt(meas(metric_func(x)))
  inv_metric_func(x) = inv(metric_func(x))

  metric = CellField(metric_func,Ω)
  sq_meas = CellField(sq_meas_func,Ω)
  inv_metric = CellField(inv_metric_func,Ω)
  MetricInfo(metric,sq_meas,inv_metric)
end



function surface_gradient(a::CellField,m::MetricInfo)
  m.inv_metric⋅ gradient(a)
end

function surface_divergence(v::CellField,m::MetricInfo)
  f = m.sq_meas*v
  1/m.sq_meas * divergence(f)
end

function surface_laplacian(f::CellField,m::MetricInfo)
  surface_divergence(surface_gradient(f,m),m)
end

m = MetricInfo(metric_func,Ω)

a = CellField(x->x[1],Ω)
b = CellField(x->VectorValue(x[1],x[2]),Ω)

surf_grad = surface_gradient(a,m)
surf_grad(pts)

surf_div = surface_divergence(b,m)
surf_div(pts)

surf_lap = surface_laplacian(a,m)
surf_lap(pts)

evaluate(surface_gradient(a,m),pts)






# XYZ = get_ambient_cell_coordinates(manifold_grid)
# panel1_3D = lazy_map(PanelMap(),XYZ, panel_ids)
# panel1_2D = lazy_map(BumpMap(), panel1_3D)

# αβ = get_cell_coordinates(manifold_grid)

# u_parametric(panel1_2D)
# map(x->u_parametric(x),panel1_2D)

# struct Ambient2ParametricFunctionMap
# end

# function Gridap.Arrays.return_cache(f::PanelRotationMap,cellx::AbstractArray{<:VectorValue})
#   A = f.mats
#   x = first(cellx)
#   T = typeof(A⋅x)
#   y = similar(cellx,T)
#   return y # CachedArray(y)
# end

# function Gridap.Arrays.evaluate!(cache,f::PanelRotationMap,cellx::AbstractArray{<:VectorValue})
#   # setsize!(cache,size(cellx))
#   y = cache
#   A = f.mats
#   map!(x -> A⋅x, y, cellx)
#   return y
# end
