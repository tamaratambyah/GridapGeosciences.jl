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
include("structs.jl")
include("refinement.jl")
include("transfer_FEmap.jl")

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
ref_cube_model = Gridap.Adaptivity.refine(cube_model)
function map_cube_to_sphere(XYZ)
  R = 1
  x,y,z = XYZ
  x_sphere = x*sqrt(1.0-y^2/2-z^2/2+y^2*z^2/(3.0))
  y_sphere = y*sqrt(1.0-z^2/2-x^2/2+x^2*z^2/(3.0))
  z_sphere = z*sqrt(1.0-x^2/2-y^2/2+x^2*y^2/(3.0))
  # R*Point(x_sphere,y_sphere,z_sphere)
  Point(2*x,4*y,z)
end

### Refine using polynomial grid
CSgrid = CubedSphereGrid(cube_grid,map_cube_to_sphere)
cmaps = collect( get_cell_map(CSgrid) )
evaluate(cmaps[1],Point(0,1))

CSmodel = CubedSphereDiscreteModel(CSgrid,topo,face_labels)
writevtk(CSmodel,datadir("CubedSphere")*"/CSmodel",append=false)

CSmodel_refined = Gridap.Adaptivity.refine(CSmodel)
writevtk(CSmodel_refined,datadir("CubedSphere")*"/CSmodel_refined",append=false)
cmaps = collect( get_cell_map(get_grid(CSmodel_refined)) )
evaluate(cmaps[1],Point(0,1))

CSmodel_refined2 = Gridap.Adaptivity.refine(CSmodel_refined)
writevtk(CSmodel_refined2,datadir("CubedSphere")*"/CSmodel_refined2",append=false)
cmaps = collect( get_cell_map(get_grid(CSmodel_refined2)) )
evaluate(cmaps[1],Point(0,1))


# interpolate map into vector FE spaces
CSgrid = CubedSphereGrid(cube_grid,map_cube_to_sphere,1)
cmaps = collect( get_cell_map(CSgrid))
evaluate(cmaps[1],Point(0,1))
CSmodel = CubedSphereDiscreteModel(CSgrid,topo,face_labels)

# writevtk(CSmodel,datadir("CubedSphere")*"/CSmodel",append=false)

CSmodel_refined = Gridap.Adaptivity.refine(CSmodel)
# writevtk(CSmodel_refined,datadir("CubedSphere")*"/CSmodel_refined",append=false)
