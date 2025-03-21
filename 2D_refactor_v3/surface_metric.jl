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

metric(x) = TensorValue{2,2}([1.0 0.0
                        0.0 1.0])
inv_metric(x) = TensorValue{2,2}([1.0 0.0
                            0.0 1.0])
sqrt_det(x) = 1.0

function u_parametric(αβ)
  α,β = αβ
  4.0.*α + β.^2
end


_model = cube_model_3D
manifold_grid = CubeGrid(_model)
panel_ids = get_panel_ids(manifold_grid)


Ω = BodyFittedTriangulation(_model, manifold_grid, IdentityVector(num_cells(manifold_grid)))
pts = get_cell_points(Ω)

a = CellField(u_parametric,Ω)
evaluate(a,pts)

"""
surface_gradient -- cellfield input
"""




function surface_gradient(a::CellField,ginv::Function)
  _ginv = CellField(ginv,get_triangulation(a))
  surface_gradient(a,_ginv)
end

function surface_gradient(a::CellField,ginv::CellField)
  ginv⋅( gradient(a) )
end

evaluate(surface_gradient(a,_inv_metric),pts)


;
# function surface_divergence(u::Function,g::Function)
#   function _surface_divergence(θϕ)
#      _g=g(θϕ)
#      function f(θϕ)
#         sqrt(det(_g))*(u(θϕ))
#      end
#      1.0/sqrt(det(_g))*(∇⋅(f))(θϕ)
#   end
# end


# function surface_laplacian(f::Function,g::Function)
#   function _surface_laplacian(θϕ)
#     surface_divergence(surface_gradient(f,g),g)(θϕ)
#   end
# end

XYZ = get_ambient_cell_coordinates(manifold_grid)
panel1_3D = lazy_map(PanelMap(),XYZ, panel_ids)
panel1_2D = lazy_map(BumpMap(), panel1_3D)

αβ = get_cell_coordinates(manifold_grid)

u_parametric(panel1_2D)
map(x->u_parametric(x),panel1_2D)

struct Ambient2ParametricFunctionMap
end

function Gridap.Arrays.return_cache(f::PanelRotationMap,cellx::AbstractArray{<:VectorValue})
  A = f.mats
  x = first(cellx)
  T = typeof(A⋅x)
  y = similar(cellx,T)
  return y # CachedArray(y)
end

function Gridap.Arrays.evaluate!(cache,f::PanelRotationMap,cellx::AbstractArray{<:VectorValue})
  # setsize!(cache,size(cellx))
  y = cache
  A = f.mats
  map!(x -> A⋅x, y, cellx)
  return y
end




# function surface_gradient(f::Function,g::Function,θϕ::SVector)
#   println("function surface grad svector")

#   surface_gradient(f,g,Point(θϕ))
# end
