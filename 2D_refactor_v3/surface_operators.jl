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
include("metric/surface_metric.jl")

parametric_model = UnstructuredDiscreteModel(CartesianDiscreteModel((-π,π,-π/2,π/2),(8,8)))
Ω = Triangulation(parametric_model)
pts = get_cell_points(Ω)


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


## auto diff
function surface_gradient(f::Number,m::MetricInfo)
  function grad_f(x::Point)
    zero(return_type(outer,x,f))
  end
end


function surface_gradient(f::Function,m::MetricInfo)
  function _gradient(x)
    m.inv_metric(x) ⋅ gradient(f,x)
  end
end

function surface_divergence(f::Function,m::MetricInfo)
  function _divergence(x)
    _f(y) =  m.sq_meas_func(y)*f(y)
    1/m.sq_meas(x) * divergence(_f,x)
  end
end


### auto diff for laplacian not working ...


# # tr(ForwardDiff.jacobian(
# x = Point(0.0,0.0)
# get_array(x)

# using ForwardDiff
# using StaticArrays
# f(x) = x[1] + x[2]


# grad(y) = ForwardDiff.gradient(f,y)

# _f(y) = convert(SMatrix{2,2,Float64},m.inv_metric_func(y) ) * grad(y)
# _f(get_array(x))

# tr(ForwardDiff.jacobian(y->_f(y), get_array(x) ))

# convert(SMatrix{2,2,Float64},m.inv_metric_func(get_array(x)) ) * grad(get_array(x))

# tr(ForwardDiff.jacobian(y->( m.inv_metric_func(y) ⋅ ForwardDiff.gradient(f,y) ), get_array(x) ))

# m.inv_metric_func(x)




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
