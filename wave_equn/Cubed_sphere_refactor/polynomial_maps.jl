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

dir = datadir("CubedSphereRefactor/Polynomial_maps")
!isdir(dir) && mkdir(dir)

cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()

### Analytical map
CSgrid_analytical_H = CubedSphereGrid(cube_grid,map_cube_to_sphere)
CSmodel_analytical_H = CubedSphereDiscreteModel(CSgrid_analytical_H,topo,face_labels)
_CSmodel_analytical_h = Gridap.Adaptivity.refine(CSmodel_analytical_H)
CSmodel_analytical_h = Gridap.Adaptivity.refine(_CSmodel_analytical_h)
writevtk(CSmodel_analytical_h,dir*"/CSmodel_analytical_h",append=false)


### Polynomial map -- without transfer
order = 1
CSgrid_poly_H = CubedSphereGrid(cube_grid,map_cube_to_sphere,order)
CSmodel_poly_H = CubedSphereDiscreteModel(CSgrid_poly_H,topo,face_labels)
_CSmodel_poly_h = Gridap.Adaptivity.refine(CSmodel_poly_H)
CSmodel_poly_h = Gridap.Adaptivity.refine(_CSmodel_poly_h)
writevtk(CSmodel_poly_h,dir*"/CSmodel_poly_h",append=false)


### Polynomial map -- with transfer
CSgrid_polyTrans_H = CubedSphereGrid(cube_grid,map_cube_to_sphere,order+1;transfer=true)
CSmodel_polyTrans_H = CubedSphereDiscreteModel(CSgrid_polyTrans_H,topo,face_labels)
_CSmodel_polyTrans_h = Gridap.Adaptivity.refine(CSmodel_polyTrans_H)
CSmodel_polyTrans_h = Gridap.Adaptivity.refine(_CSmodel_polyTrans_h)
writevtk(CSmodel_polyTrans_h,dir*"/CSmodel_polyTrans_h",append=false)


using GridapGeosciences
alberto = GridapGeosciences.CubedSphereDiscreteModel(4,order)
writevtk(alberto,dir*"/alberto",append=false)
