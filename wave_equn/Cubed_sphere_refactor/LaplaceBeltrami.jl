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
using GridapGeosciences

include("cube_surface_1_cell_per_panel.jl")
include("structs.jl")
include("refinement.jl")
include("maps.jl")
include("normals.jl")

dir = datadir("CubedSphereRefactor/Laplace")
!isdir(dir) && mkdir(dir)

function _divergence_unit_sphere(v)
  function tmp(θϕ)
     Gθϕ=GridapGeosciences.G_unit_sphere(θϕ)
     Jθϕ=GridapGeosciences.J_unit_sphere(θϕ)
     function f(θϕ)
        sqrt(det(Gθϕ))*(transpose(Jθϕ)⋅(v(θϕ))) # convert to ambient
     end
       1.0/sqrt(det(Gθϕ))*(∇⋅(f))(θϕ)
  end
end

function _laplacian_unit_sphere(v)
  function tmp(θϕ)
    _divergence_unit_sphere(gradient_unit_sphere(v))(θϕ)
  end
end

function u(x)
  x[1]*x[2]*x[3]
  # sin(π*x[1])*cos(π*x[2])*exp(x[3])
end
function uθϕ(θϕ)
  u(θϕ2xyz(θϕ))
end
function f(x)
  -(_laplacian_unit_sphere(uθϕ)(xyz2θϕ(x)))
end
function Gridap.∇(::typeof(u))
  gradient_unit_sphere(uθϕ)
end


cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()
CSgrid_analytical_H = CubedSphereGrid(cube_grid,map_cube_to_sphere)
CSmodelH = CubedSphereDiscreteModel(CSgrid_analytical_H,topo,face_labels)
_CSmodel = Gridap.Adaptivity.refine(CSmodelH)
CSmodel = Gridap.Adaptivity.refine(_CSmodel)

l2(w,dΩ) = sum( ∫( w⊙w )dΩ  )

model = CSmodel
order = 3
reffe = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,reffe;conformity=:H1)
U = TrialFESpace(V)
degree = 2*order
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

a(u,v) = ∫( ∇(v)⋅∇(u) )*dΩ
b(v) = ∫( v*f )*dΩ
op = AffineFEOperator(a,b,U,V)
uh = solve(LUSolver(),op)
e = u-uh
l2(e,dΩ)

writevtk(Ω,dir*"/manufactured_sol",cellfields=["u"=>u,"uh"=>uh,"e"=>e],append=false)
