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
include("cube_sphere_1_cell_per_panel.jl")
include("structs.jl")
include("refinement.jl")
include("maps.jl")

dir = datadir("CubedSphereRefactor")
!isdir(dir) && mkdir(dir)


cube_grid,topo,face_labels = cube_sphere_1_cell_per_panel()
cube_model = UnstructuredDiscreteModel(cube_grid,topo,face_labels)
ref_cube_model = Gridap.Adaptivity.refine(cube_model)
cube_nodes = get_node_coordinates(cube_grid)

### Analytical map
CSgrid = CubedSphereGrid(cube_grid,map_cube_to_sphere)
CSmodel = CubedSphereDiscreteModel(CSgrid,topo,face_labels)
writevtk(CSmodel,dir*"/CSmodel",append=false)

# test node 4 = (1,1,1) is mapped properly in cells 1,2,3
node4 = cube_nodes[4] # node 4 = (1,1,1)
mapped_node4 = map_cube_to_sphere(node4)
cmaps = collect( get_cell_map(CSmodel) )
@test evaluate(cmaps[1],Point(1,1)) == mapped_node4
@test evaluate(cmaps[2],Point(1,0)) == mapped_node4
evaluate(cmaps[3],Point(0,0)) == mapped_node4

# apply refinement
CSmodel_refined1 = Gridap.Adaptivity.refine(CSmodel)
writevtk(CSmodel_refined1,dir*"/CSmodel_refined1",append=false)

CSmodel_refined2 = Gridap.Adaptivity.refine(CSmodel_refined1)
writevtk(CSmodel_refined2,dir*"/CSmodel_refined2",append=false)


### Polynomial map
CSgrid_poly = CubedSphereGrid(cube_grid,map_cube_to_sphere,1)
CSmodel_poly = CubedSphereDiscreteModel(CSgrid_poly,topo,face_labels)
writevtk(CSmodel_poly,dir*"/CSmodel_poly",append=false)

# test node 4 = (1,1,1) is mapped properly in cells 1,2,3 for mapping of order 1
cmaps_poly = collect( get_cell_map(CSmodel_poly) )
@test evaluate(cmaps_poly[1],Point(1,1)) == mapped_node4
@test evaluate(cmaps_poly[2],Point(1,0)) == mapped_node4
@test evaluate(cmaps_poly[3],Point(0,0)) == mapped_node4

# apply order 2 mapping
CSgrid_poly2 = CubedSphereGrid(cube_grid,map_cube_to_sphere,2)
CSmodel_poly2 = CubedSphereDiscreteModel(CSgrid_poly2,topo,face_labels)

# apply refinement
CSmodel_poly2_refined1 = Gridap.Adaptivity.refine(CSmodel_poly2)
writevtk(CSmodel_poly2_refined1,dir*"/CSmodel_poly_refined1",append=false)

CSmodel_poly_refined2 = Gridap.Adaptivity.refine(CSmodel_poly2_refined1)
writevtk(CSmodel_poly2_refined2,dir*"/CSmodel_poly_refined2",append=false)


using GridapGeosciences
alberto = GridapGeosciences.CubedSphereDiscreteModel(4,2; radius=1)
writevtk(alberto,dir*"/alberto",append=false)
