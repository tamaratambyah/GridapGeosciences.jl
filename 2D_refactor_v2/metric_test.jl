using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Adaptivity
using Test
using LinearAlgebra
using FillArrays
using BenchmarkTools
include("cube_topo/cube_surface_1_cell_per_panel.jl")

metric = [1 0
          0 1]

g = map(x->TensorValue(metric),1:6)

_g(x) = TensorValue(metric)





"""
surface_gradient -- function input
"""
function surface_gradient(f::Function,g::Function)
  println("function surface grad")
  function _surface_gradient(θϕ)
    surface_gradient(f,g,θϕ)
  end
end

function surface_gradient(f::Function,g::Function,θϕ::Point)
  println("function surface grad point")

  _g=g(θϕ)
  inv(_g)⋅( gradient(f,θϕ) )
end

function surface_gradient(f::Function,g::Function,θϕ::SVector)
  println("function surface grad svector")

  surface_gradient(f,g,Point(θϕ))
end

"""
surface_gradient -- cellfield input
"""
function surface_gradient(a::CellField,g::Function)
  println("cellfield function surface grad")
  _g = CellField(g,get_triangulation(a))
  surface_gradient(a,_g)
end

function surface_gradient(a::CellField,g::CellField)
  println("cellfield cellfield surface grad")

  println(inv(g))
  ( gradient(a) )
  # inv(g)⋅( gradient(a) )
end


function surface_divergence(u::Function,g::Function)
  function _surface_divergence(θϕ)
     _g=g(θϕ)
     function f(θϕ)
        sqrt(det(_g))*(u(θϕ))
     end
     1.0/sqrt(det(_g))*(∇⋅(f))(θϕ)
  end
end


function surface_laplacian(f::Function,g::Function)
  function _surface_laplacian(θϕ)
    surface_divergence(surface_gradient(f,g),g)(θϕ)
  end
end
