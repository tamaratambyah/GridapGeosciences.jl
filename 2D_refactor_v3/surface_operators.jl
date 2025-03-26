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
include("surface_metric_and_op/operators.jl")
include("surface_metric_and_op/quadrature.jl")
include("surface_metric_and_op/cubedsphere_metric.jl")


r = sqrt(3.0)
parametric_model = UnstructuredDiscreteModel(CartesianDiscreteModel((-π/4,π/4,-π/4,π/4),(8,8)))
Ω = Triangulation(parametric_model)
pts = get_cell_points(Ω)

order = 4
m = MetricInfo(metric_func,Ω)
dΩg = Measure(m,Ω,order)

a = CellField(x->x[1],Ω)
b = CellField(x->VectorValue(x[1],x[2]),Ω)

surf_grad = surface_gradient(a,m)
surf_grad(pts)

surf_div = surface_divergence(b,m)
surf_div(pts)

surf_lap = surface_laplacian(a,m)
surf_lap(pts)

evaluate(surface_gradient(a,m),pts)

#### try combining with surface operators

f(x) = x[1] + x[2]
h(x) = VectorValue(2.0*x[1], 3.0)
sum( ∫( surface_gradient(f,m)  )dΩg )

sum( ∫( surface_divergence(h,m) )dΩg )


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
