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


cube_grid,topo,face_labels = cube_sphere_1_cell_per_panel()

# mapp(x) = Point(2*x[1],4*x[2],-1*x[3])
mapp(x) = map_cube_to_sphere(x)
order = 3

# coarse model, space, and FE map
cube_modelH = UnstructuredDiscreteModel(cube_grid)
VH = FESpace(cube_modelH,
            ReferenceFE(lagrangian,VectorValue{3,Float64},order),
            conformity=:H1)
mapH = interpolate(mapp,VH)

# refined model, space, and FE map
cube_modelh = Gridap.Adaptivity.refine(cube_modelH)
Vh = FESpace(cube_modelh,
            ReferenceFE(lagrangian,VectorValue{3,Float64},order),
            conformity=:H1)
Πmaph = interpolate(mapp,Vh)

# transfer coarse → fine map
mapHh, order, _Vh = transfer_FE_map(cube_modelh,mapH)

# test error between projection of map into refined FE space and transfer of map
l2(w,dΩ) = sum( ∫( w⊙w )dΩ  )
Ωh = Triangulation(cube_modelh)
dΩh = Measure(Ωh,2*1+1)
error = l2(Πmaph - mapHh,dΩh)
@test error < 1e-10
