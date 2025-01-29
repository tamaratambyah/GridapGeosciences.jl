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

function map_parameteric_sphere(XYZ)
  sphereXYZ = map_cube_to_sphere(XYZ)
  GridapGeosciences.xyz2θϕ(sphereXYZ)
end

x = Point(1,1,1)
X = map_cube_to_sphere(x)
θϕ = GridapGeosciences.xyz2θϕ(X)

(GridapGeosciences.xyz2θϕ∘map_cube_to_sphere)(x)
map_parameteric_sphere(x)





dir = datadir("CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()
cube_model = UnstructuredDiscreteModel(cube_grid,topo,face_labels)
ref_cube_model = Gridap.Adaptivity.refine(cube_model)
cube_nodes = get_node_coordinates(cube_grid)

### Analytical map
CSgrid = CubedSphereGrid(cube_grid,map_parameteric_sphere)
CSmodel = CubedSphereDiscreteModel(CSgrid,topo,face_labels)
writevtk(CSmodel,dir*"/CSmodel",append=false)

# test node 4 = (1,1,1) is mapped properly in cells 1,2,3
node4 = cube_nodes[4] # node 4 = (1,1,1)
mapped_node4 = map_parameteric_sphere(node4)
cmaps = collect( get_cell_map(CSmodel) )
@test evaluate(cmaps[1],Point(1,1)) == mapped_node4
@test evaluate(cmaps[2],Point(1,0)) == mapped_node4
@test evaluate(cmaps[3],Point(0,0))== mapped_node4

# apply refinement
CSmodel_refined1 = Gridap.Adaptivity.refine(CSmodel)
writevtk(CSmodel_refined1,dir*"/CSmodel_refined1",append=false)

CSmodel_refined2 = Gridap.Adaptivity.refine(CSmodel_refined1)
writevtk(CSmodel_refined2,dir*"/CSmodel_refined2",append=false)

abstract type CubedSphereDiscreteModel{Dc,Dp,_Dp} <: Grid{Dc,Dp} end

struct CSDiscreteModel{Dc,Dp} <: CubedSphereDiscreteModel{Dc,Dp,_Dp}
  data
end

struct AdaptedCSDiscreteModel{Dc,Dp} <: CubedSphereDiscreteModel{Dc,Dp,_Dp}
  data
end
