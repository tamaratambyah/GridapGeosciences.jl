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
CSmodel = CubedSphereDiscreteModel(CSgrid_analytical_H,topo,face_labels)

# plot the normals
normals = get_outward_normal_vector(CSmodel)
writevtk(Triangulation(CSmodel),dir*"/normals",cellfields=["n"=>normals],append=false)


### Low level construction of normals
cell_Jt,cell_x = get_tangent_space_basis(CSmodel)

# check the normal of the first cell is (1,0,0,)
_v = map(evaluate,cell_Jt[1],cell_x)
v = _v[1][1]

# cross product of basis vectors
n1 = v[1,2]*v[2,3] - v[1,3]*v[2,2]
n2 = v[1,3]*v[2,1] - v[1,1]*v[2,3]
n3 = v[1,1]*v[2,2] - v[1,2]*v[2,1]
n = VectorValue(n1,n2,n3)
@test n/norm(n) == Point(1,0,0)

# print all the normals
cell_normal = lazy_map( Operation(_unit_outward_normal), cell_Jt)
map(cell_normal) do n
  println( map(evaluate,n,cell_x) )
end
nx = map(evaluate,cell_normal[1],cell_x)
