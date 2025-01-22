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

dir = datadir("CubedSphereRefactor/Normals")
!isdir(dir) && mkdir(dir)


### Analytical map
cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()
CSgrid_analytical_H = CubedSphereGrid(cube_grid,map_cube_to_sphere)
_CSmodel = CubedSphereDiscreteModel(CSgrid_analytical_H,topo,face_labels)
CSmodel = Gridap.Adaptivity.refine(_CSmodel)

Ω = Triangulation(CSmodel)
Λ = SkeletonTriangulation(CSmodel)
writevtk(Λ,dir*"/skeleton",append=false)

n_Λ = get_normal_vector(Λ)

u = CellField(Point(1,1,1),Ω)


writevtk(Λ,dir*"/skeleton",
        cellfields=["plus"=>n_Λ.plus,
                    "minus"=>n_Λ.minus,
                    "u1"=>u⋅(n_Λ.plus),
                    "u2"=>u⋅(n_Λ.minus)],append=false)


topo = get_grid_topology(CSmodel)
D = num_cell_dims(CSmodel)
face_to_mask = collect(Bool, .!get_isboundary_face(topo,D-1))
SkeletonTriangulation(model,face_to_mask)
