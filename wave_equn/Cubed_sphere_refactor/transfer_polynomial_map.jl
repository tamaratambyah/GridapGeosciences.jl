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


dir = datadir("CubedSphereRefactor/Transfer_polynomial_maps")
!isdir(dir) && mkdir(dir)

cube_grid,topo,face_labels = cube_surface_1_cell_per_panel()
cube_model0 = UnstructuredDiscreteModel(cube_grid)
cube_modelH = Gridap.Adaptivity.refine(cube_model0)
cube_modelh = Gridap.Adaptivity.refine(cube_modelH)

writevtk(cube_modelH,dir*"/cube_modelH",append=false)
writevtk(cube_modelh,dir*"/cube_modelh",append=false)

l2(w,dΩ) = sum( ∫( w⊙w )dΩ  )
mapp(x) = map_cube_to_sphere(x)
# mapp(x) = map_cube_to_cube(x)

for order = [1,2,3,4]

  # coarse model, space, and FE map
  VH = FESpace(cube_modelH,
              ReferenceFE(lagrangian,VectorValue{3,Float64},order),
              conformity=:H1)
  mapH = interpolate(mapp,VH)


  # refined model, space, and FE map
  Vh = FESpace(cube_modelh,
              ReferenceFE(lagrangian,VectorValue{3,Float64},order),
              conformity=:H1)
  Πmaph = interpolate(mapp,Vh)

  # transfer coarse → fine map
  mapHh, _Vh = transfer_FE_map(cube_modelh,mapH,order)

  # test error between projection of map into refined FE space and transfer of map
  Ωh = Triangulation(cube_modelh)
  dΩh = Measure(Ωh,4*order)
  error = l2(Πmaph - mapHh,dΩh)
  println("order = ", order, "; error = ",error)
end
