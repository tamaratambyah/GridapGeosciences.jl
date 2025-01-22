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
include("refinement.jl")
include("maps.jl")
include("normals.jl")

dir = datadir("CubedSphereRefactor/Laplace")
!isdir(dir) && mkdir(dir)

# cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()
# CSgrid_analytical_H = CubedSphereGrid(cube_grid,map_cube_to_sphere)
# CSmodel = CubedSphereDiscreteModel(CSgrid_analytical_H,topo,face_labels)


u(x) = x[1]*(1-x[1]) + x[2]*(1-x[2])
f(x) = -laplacian(u)(x) #+ 2
l2(w,dΩ) = sum( ∫( w⊙w )dΩ  )

model = CartesianDiscreteModel((0,1,0,1),(10,10))
order = 2
reffe = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,reffe;conformity=:H1,dirichlet_tags="boundary")
U = TrialFESpace(V,u)
degree = 2*order
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

a(u,v) = ∫( ∇(v)⋅∇(u) )*dΩ
b(v) = ∫( v*f )*dΩ
op = AffineFEOperator(a,b,U,V)
uh = solve(LUSolver(),op)
e = u-uh
l2(e,dΩ)
