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
using StaticArrays: SVector

include("cube_surface_1_cell_per_panel.jl")
include("structs.jl")
include("refinement.jl")
include("maps.jl")


# metric(θϕ)
u(θϕ) = 3*θϕ[1] + 2*θϕ[2]^2
w(θϕ) = VectorValue( 3*θϕ[1],2*θϕ[2])

# θϕ = VectorValue(π,π/4)
# evaluate(surface_gradient(u,metric),θϕ)
# evaluate(surface_divergence(w),θϕ)
# evaluate(surface_laplacian(u),θϕ)

function metric(θϕ)
  θ,ϕ = θϕ
  # TensorValue{2,2}( (cos(ϕ))^2, 0,
  #                     0,        1)

  TensorValue{2,2}( 1, 0, 0, 1)
end

θs = VectorValue(0,2π)
ϕs = VectorValue(-π/2,π/2)
θϕ = VectorValue(θs[2],ϕs[1])

g = metric(θϕ)
invg = inv(g)

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

f(θϕ) = surface_laplacian(u,metric)(θϕ)


dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()

CSgrid = CubedSphereGrid(cube_grid,map_cube_to_latlon)

CSmodel = ManifoldDiscreteModel(CSgrid,topo,face_labels)
CSmodelh = Gridap.Adaptivity.refine(CSmodel)
model = Gridap.Adaptivity.refine(CSmodelh)

writevtk(model,dir*"/model",append=false)


p = 1
degree = 2*p+1
V = FESpace(model,ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

_metric = CellField(metric,Ω)
invg = evaluate(inv(_metric),get_cell_points(Ω))

pt = (get_node_coordinates(Ω))


_f = CellField(f,Ω)
a(u,v) = ∫( surface_gradient(u,metric)⋅surface_gradient(v,metric) )dΩ
b(v)   = ∫( 0.0*v )dΩ

res(u,v) = a(u,v) - b(v)
jac(u,du,v) = a(du,v)

op = FEOperator(res,jac,U,V)
solver = NLSolver(LUSolver())
uh = solve(solver,op)

e=u-uh
∇e = ∇(uh)-∇(u)∘xyz2θϕ

# H1-norm
# sqrt(sum(∫( e*e + (∇e)⋅(∇e) )dΩ))
sqrt(sum(∫( (∇e)⋅(∇e) )dΩ))#,uh
