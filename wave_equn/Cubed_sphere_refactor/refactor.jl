using DrWatson
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Test
using LinearAlgebra
using FillArrays


nodes = [
Point(1.0, -1.0, -1.0)
Point(1.0, 1.0, -1.0)
Point(1.0, -1.0, 1.0)
Point(1.0, 1.0, 1.0)
Point(-1.0, -1.0, 1.0)
Point(-1.0, 1.0, 1.0)
Point(-1.0, 1.0, -1.0)
Point(-1.0, -1.0, -1.0)
]

function generate_ptr(n)
  nvertices = 4
  ptr  = Vector{Int}(undef,n+1)
  ptr[1]=1
  for i=1:n
    ptr[i+1]=ptr[i]+nvertices
  end
  ptr
end

data = [ 1,2,3,4, 3,4,5,6, 4,2,6,7, 6,7,5,8, 7,2,8,1, 8,1,5,3  ]
ptr = generate_ptr(6)
cell_node_ids = Table(data,ptr)


polytopes = fill(QUAD,6)
cell_type = fill(1,6)
reffes = LagrangianRefFE(Float64,QUAD,1)
cell_reffes=[reffes]

topo = UnstructuredGridTopology(nodes,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
face_labels = FaceLabeling(topo)

cube_grid = Gridap.Geometry.UnstructuredGrid(nodes,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

cube_model = UnstructuredDiscreteModel(cube_grid,topo,face_labels)
# refined_cube_model = Gridap.Adaptivity.refine(cube_model)

# refined_cube_model = Gridap.Adaptivity.refine(cube_model)
# writevtk(refine_model,datadir("CubedSphere")*"/refine_model",append=false)


CSgrid = CubedSphereGrid(cube_grid,map_cube_to_sphere)
CSmodel = DiscreteModel(CSgrid,topo,face_labels)
_CSmodel = UnstructuredDiscreteModel(CSmodel)

CSmodel_refined = Gridap.Adaptivity.refine(_CSmodel)
writevtk(CSmodel_refined,datadir("CubedSphere")*"/CSmodel_refined1",append=false)

## Refined
refine_cube_grid = get_grid(refined_cube_model)
refine_topo = get_grid_topology(refined_cube_model)
refine_face_labels = get_face_labeling(refined_cube_model)

CSgrid = CubedSphereGrid(refine_cube_grid,map_cube_to_sphere)
cmaps = collect( get_cell_map(CSgrid) )
evaluate(cmaps[1],Point(0,1))




function map_cube_to_sphere(XYZ)
  x,y,z = XYZ
  x_sphere = x*sqrt(1.0-y^2/2-z^2/2+y^2*z^2/(3.0))
  y_sphere = y*sqrt(1.0-z^2/2-x^2/2+x^2*z^2/(3.0))
  z_sphere = z*sqrt(1.0-x^2/2-y^2/2+x^2*y^2/(3.0))
  Point(x_sphere,y_sphere,z_sphere)
  # Point(2*x,4*y,z)
end



# Dc = num_cell_dims(topo) = dimension of model = 2
# Dp = num_point_dims(grid)  = dimension of points = 3
# Tp = VectorValue{3,Float64}
# B =

struct CubedSphereGrid <: Grid{2,3}
  cube_grid
  sphere_grid
  sphere_cell_map
end
# Dc = 2
# Dp = 3
# Tp = VectorValue{3,Float64}
# O = typeof( Gridap.Geometry.NonOriented() )
# Tn = typeof( nothing )
# struct _CubedSphereGrid{Dc,Dp,Tp,O,Tn} <: UnstructuredGrid{Dc,Dp,Tp,O,Tn}
#   cube_grid
#   sphere_grid
#   sphere_cell_map
# end

function CubedSphereGrid(cube_grid,map_cube_to_sphere::Function)

  # map the nodes from the cube to the sphere
  cube_nodes = get_node_coordinates(cube_grid)
  sphere_nodes = map( (XYZ)-> map_cube_to_sphere(XYZ), cube_nodes  ) #lazy_map( map_cube_to_sphere, cube_nodes  )

  sphere_grid = Gridap.Geometry.UnstructuredGrid(sphere_nodes,
                            get_cell_node_ids(cube_grid),
                            get_reffes(cube_grid),
                            get_cell_type(cube_grid),
                            Gridap.Geometry.NonOriented())


  # create map from 2D reference FE to the sphere
  geo_cell_map = Fill( Gridap.Fields.GenericField( map_cube_to_sphere ), num_cells(cube_grid))
  cube_cell_map = get_cell_map(cube_grid)
  sphere_cell_map = lazy_map(∘,geo_cell_map,cube_cell_map)

  CubedSphereGrid(cube_grid,sphere_grid,sphere_cell_map)

end

Gridap.Geometry.OrientationStyle(grid::CubedSphereGrid) = Gridap.Geometry.NonOriented()

function Gridap.Geometry.get_reffes(grid::CubedSphereGrid)
  Gridap.Geometry.get_reffes(grid.sphere_grid)
end

function Gridap.Geometry.get_cell_type(grid::CubedSphereGrid)
  Gridap.Geometry.get_cell_type(grid.sphere_grid)
end

function Gridap.Geometry.get_node_coordinates(grid::CubedSphereGrid)
  Gridap.Geometry.get_node_coordinates(grid.sphere_grid)
end

function Gridap.Geometry.get_cell_node_ids(grid::CubedSphereGrid)
  Gridap.Geometry.get_cell_node_ids(grid.sphere_grid)
end

function Gridap.Geometry.get_cell_map(grid::CubedSphereGrid)
  grid.sphere_cell_map
end

### Original
CSgrid = CubedSphereGrid(cube_grid,map_cube_to_sphere)
CSmodel = UnstructuredDiscreteModel( DiscreteModel(CSgrid,topo,face_labels) )
writevtk(CSmodel,datadir("CubedSphere")*"/CSmodel",append=false)

## Refined
refine_cube_grid = get_grid(refined_cube_model)
refine_topo = get_grid_topology(refined_cube_model)
refine_face_labels = get_face_labeling(refined_cube_model)

CSgrid = CubedSphereGrid(refine_cube_grid,map_cube_to_sphere)
cmaps = collect( get_cell_map(CSgrid) )
evaluate(cmaps[1],Point(0,1))


CSmodel_refined = UnstructuredDiscreteModel( DiscreteModel(CSgrid,refine_topo,refine_face_labels) )
writevtk(CSmodel_refined,datadir("CubedSphere")*"/CSmodel_refined",append=false)

# refined_CSmodel = Gridap.Adaptivity.refine(CSmodel)
# writevtk(refined_CSmodel,datadir("CubedSphere")*"/refined_CSmodel",append=false)
