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


dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()

CSgrid = CubedSphereGrid(cube_grid,map_cube_to_latlon)
CSmodel = ManifoldDiscreteModel(CSgrid,topo,face_labels)

CSmodelh = Gridap.Adaptivity.refine(CSmodel)
CSmodelh2 = Gridap.Adaptivity.refine(CSmodelh)

writevtk(CSmodel,dir*"/CSmodel",append=false)
writevtk(CSmodelh,dir*"/CSmodelh",append=false)
writevtk(CSmodelh2,dir*"/CSmodelh2",append=false)


cube_nodes = get_node_coordinates(cube_grid)


node4 = cube_nodes[4] # node 4 = (1,1,1)
mapped_node4 = map_cube_to_latlon(node4)
cmaps = collect(get_cell_map(CSmodel))
@test evaluate(cmaps[1],Point(1,1)) == mapped_node4
@test evaluate(cmaps[2],Point(1,0)) == mapped_node4
evaluate(cmaps[3],Point(0,0)) == mapped_node4
