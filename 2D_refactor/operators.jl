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
include("cube_surface_1_cell_per_panel.jl")
include("structs.jl")
include("maps.jl")

θϕ = VectorValue(π,π/4)

metric(θϕ)
u(θϕ) = 3*θϕ[1] + 2*θϕ[2]^2
w(θϕ) = VectorValue( 3*θϕ[1],2*θϕ[2])

evaluate(surface_gradient(u),θϕ)
evaluate(surface_divergence(w),θϕ)
evaluate(surface_laplacian(u),θϕ)

function metric(θϕ)
  θ,ϕ = θϕ
  TensorValue{2,2}( (cos(ϕ))^2, 0,
                      0,        1)
end


function surface_gradient(u)
  function _surface_gradient(θϕ)
     g=metric(θϕ)
     inv(g)⋅(∇(u))(θϕ)
  end
end


function surface_divergence(u)
  function _surface_divergence(θϕ)
     g=metric(θϕ)
     function f(θϕ)
        sqrt(det(g))*(u(θϕ))
     end
     1.0/sqrt(det(g))*(∇⋅(f))(θϕ)
  end
end


function surface_laplacian(u)
  function _surface_laplacian(θϕ)
    surface_divergence(surface_gradient(u))(θϕ)
  end
end

f(θϕ) = surface_laplacian(u)(θϕ)


dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()

CSgrid = CubedSphereGrid(cube_grid,map_cube_to_latlon)
CSmodel = ManifoldDiscreteModel(CSgrid,topo,face_labels)

Dp = num_point_dims(CSmodel)
num_point_dims(CSgrid)


V = FESpace(CSmodel,ReferenceFE(lagrangian,Float64,1); conformity=:H1)
U = TrialFESpace(V)

Ω = Triangulation(CSmodel)
dΩ = Measure(Ω,4)

a(u,v) = ∫(surface_gradient(v)⋅surface_gradient(u))dΩ
b(v)   = ∫(v*0.0)dΩ

op = AffineFEOperator(a,b,U,V)
fels = LinearFESolver(ls)
uh = solve(fels,op)

e=u-uh
∇e = ∇(uh)-∇(u)∘xyz2θϕ

# H1-norm
# sqrt(sum(∫( e*e + (∇e)⋅(∇e) )dΩ))
sqrt(sum(∫( (∇e)⋅(∇e) )dΩ))#,uh
