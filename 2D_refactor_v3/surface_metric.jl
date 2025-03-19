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

metric() = TensorValue{2,2}([1.0 0.0
                        0.0 1.0])
inv_metrix() = TensorValue{2,2}([1.0 0.0
                            0.0 1.0])

_metric(x) = TensorValue{2,2}([x[1] 0.0
                          0.0 1.0])
_inv_metric(x) = TensorValue{2,2}([0.0 x[1]
                              0.0 1.0])

_model = cube_model_3D
manifold_grid = CubeGrid(_model)
panel_ids = get_panel_ids(manifold_grid)

XYZ = get_ambient_cell_coordinates(manifold_grid)
panel1_3D = lazy_map(PanelMap(),XYZ, panel_ids)
panel1_2D = lazy_map(BumpMap(), panel1_3D)

map(x->u_parametric(x),panel1_2D[1])

function u_parametric(αβ)
  α,β = αβ
  4.0*α + β^2
end

trian = GenericTriangulation(manifold_grid)
pts = get_cell_points(trian)

a = CellField(u_parametric,trian)
evaluate(a,pts)

"""
surface_gradient -- cellfield input
"""
function surface_gradient(a::CellField,ginv::Function)
  println("cellfield function surface grad")
  _ginv = CellField(ginv,get_triangulation(a))
  surface_gradient(a,_ginv)
end

function surface_gradient(a::CellField,ginv::CellField)
  println("cellfield cellfield surface grad")
  ginv⋅( gradient(a) )
end

evaluate(surface_gradient(a,_inv_metric),pts)


G = CellField(_metric,get_triangulation(a))
evaluate(G,pts)

ginv = CellField(_inv_metric,get_triangulation(a))
evaluate(ginv,pts)

gradu = gradient(a)
evaluate(gradu,pts)

sufgrad = ginv.⋅gradu
evaluate(sufgrad,pts)







# """
# surface_gradient -- function input
# """
# function surface_gradient(f::Function,g::Function)
#   println("function surface grad")
#   function _surface_gradient(θϕ)
#     surface_gradient(f,g,θϕ)
#   end
# end

# function surface_gradient(f::Function,g::Function,θϕ::Point)
#   println("function surface grad point")

#   _g=g(θϕ)
#   inv(_g)⋅( gradient(f,θϕ) )
# end

# function surface_gradient(f::Function,g::Function,θϕ::SVector)
#   println("function surface grad svector")

#   surface_gradient(f,g,Point(θϕ))
# end



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
