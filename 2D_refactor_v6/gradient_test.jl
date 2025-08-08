using DrWatson
using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity
using LinearAlgebra

include("panel_ids_from_refinement.jl")

dir = datadir("2D_CubedSphereRefactor")
!isdir(dir) && mkdir(dir)

function generate_ptr(n)
  nvertices = 4
  ptr  = Vector{Int}(undef,n+1)
  ptr[1]=1
  for i=1:n
    ptr[i+1]=ptr[i]+nvertices
  end
  ptr
end


npanels = 2

nodes = [
  Point(0.0, 0.0)  # node 1
  Point(1.0, 0.0)   # node 2
  Point(0.0, 1.0)   # node 3
  Point(1.0, 1.0)    # node 4
  Point(2.0, 0.0)  # node 5
  Point(2.0, 1.0)   # node  6
]

## CCAM panel ordering
# data = [ 1,2,3,4, 2,5,4,6 ]
data = [ 1,2,3,4, 6,5,4,2 ]

ptr = generate_ptr(npanels)
cell_node_ids = Table(data,ptr)

polytopes = fill(QUAD,npanels)
cell_type = fill(1,npanels)
reffes = LagrangianRefFE(Float64,QUAD,1)
cell_reffes=[reffes]

topo = UnstructuredGridTopology(nodes,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
labels = FaceLabeling(topo)

grid = Gridap.Geometry.UnstructuredGrid(nodes,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

model = UnstructuredDiscreteModel(grid,topo,labels)

model = Gridap.Adaptivity.refine(model)
panel_ids = get_panel_ids(model)

u(x) = x[1]*x[2] + 2*x[1]^2
function uex(p)
  function _u(x)
    u(x)
  end
end
∇u(x) = VectorValue(x[2] + 4*x[1],x[1])

Ω = Triangulation(model)

cell_field = map(p->GenericField(uex(p)),panel_ids)
ucf =  CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

gu = gradient(ucf)

writevtk(Ω,dir*"/gradient_test",cellfields=["u"=>u,"du"=>∇u,"ucf"=>ucf,"gradu"=>gu],append=false)
