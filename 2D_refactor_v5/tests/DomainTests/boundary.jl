using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields

dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

a = 1
npanels = 6

nodes_3d = a.* [
  Point(1.0, -1.0, -1.0)  # node 1
  Point(1.0, 1.0, -1.0)   # node 2
  Point(1.0, -1.0, 1.0)   # node 3
  Point(1.0, 1.0, 1.0)    # node 4
  Point(-1.0, -1.0, 1.0)  # node 5
  Point(-1.0, 1.0, 1.0)   # node 6
  Point(-1.0, 1.0, -1.0)  # node 7
  Point(-1.0, -1.0, -1.0) # node 8
]

## CCAM panel ordering
data = [ 1,2,3,4, 3,4,5,6, 2,7,4,6, 8,5,7,6, 1,8,2,7, 1,3,8,5  ] # reorient + rotated

ptr = generate_ptr(npanels)
cell_node_ids = Table(data,ptr)

polytopes = fill(QUAD,npanels)
cell_type = fill(1,npanels)
reffes = LagrangianRefFE(Float64,QUAD,1)
cell_reffes=[reffes]

topo = UnstructuredGridTopology(nodes_3d,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())


function Gridap.Geometry.compute_isboundary_face(g::GridTopology,d::Integer)
  D = num_cell_dims(g)
  nfaces = num_faces(g,d)

  if d == D-1
    println("if")
    # Geometry._compute_isboundary_facet!(g)
    return fill(true,nfaces)
  else
    println("else")
    # Geometry._compute_isboundary_face!(g,d)
    return fill(false,nfaces)
  end
end
get_isboundary_face(topo,0)

_labels = FaceLabeling(topo)

cube_grid = Gridap.Geometry.UnstructuredGrid(nodes_3d,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

cube_model_3D = UnstructuredDiscreteModel(cube_grid,topo,_labels)

cube_model_3D = Gridap.Adaptivity.refine(cube_model_3D)


Ω = Triangulation(cube_model_3D)
Γ = BoundaryTriangulation(cube_model_3D;tags="boundary")
writevtk(Γ,dir*"/CS",append=false)
writevtk(cube_model_3D,dir*"/CS",append=false)

Measure(Γ,2)
