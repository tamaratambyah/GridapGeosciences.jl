"""
build panel 1 and 2 by hand
"""

using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../../../src/initialise.jl")

function coarse_panel1_3D(a::Real)
  npanels = 1

  nodes_3d = a.* [
    Point(1.0, -1.0, -1.0)  # node 1
    Point(1.0, 1.0, -1.0)   # node 2
    Point(1.0, -1.0, 1.0)   # node 3
    Point(1.0, 1.0, 1.0)    # node 4
    # Point(-1.0, -1.0, 1.0)  # node 5
    # Point(-1.0, 1.0, 1.0)   # node 6
  ]

  ## CCAM panel ordering
  data = [ 1,2,3,4 ] # reorient + rotated

  ptr = generate_ptr(npanels)
  cell_node_ids = Table(data,ptr)

  polytopes = fill(QUAD,npanels)
  cell_type = fill(1,npanels)
  reffes = LagrangianRefFE(Float64,QUAD,1)
  cell_reffes=[reffes]

  topo = UnstructuredGridTopology(nodes_3d,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())

  d_to_num_dfaces = [ num_faces(topo,d) for d in 0:num_dims(topo)]
  labels = FaceLabeling(d_to_num_dfaces)

  get_face_entity(labels,0) .= get_isboundary_face(topo,0) #.+ 1
  get_face_entity(labels,1) .= fill(true,length(get_isboundary_face(topo,1))) .+ 1
  get_face_entity(labels,2) .= get_isboundary_face(topo,2) #.+ 1

  add_tag!(labels,"interior",[1,])
  add_tag!(labels,"boundary",[2,])
  add_tag_from_tags!(labels,"all",["interior","boundary"])

  panel_grid = Gridap.Geometry.UnstructuredGrid(nodes_3d,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

  panel_ids = collect(1)
  return panel_grid,topo,labels,panel_ids
end

function coarse_panel2_3D(a::Real)
  npanels = 1

  nodes_3d = a.* [
    # Point(1.0, -1.0, -1.0)  # node 1
    # Point(1.0, 1.0, -1.0)   # node 2
    Point(1.0, -1.0, 1.0)   # node 3
    Point(1.0, 1.0, 1.0)    # node 4
    Point(-1.0, -1.0, 1.0)  # node 5
    Point(-1.0, 1.0, 1.0)   # node 6
  ]

  ## CCAM panel ordering
  data = [ 1,2,3,4  ] # reorient + rotated

  ptr = generate_ptr(npanels)
  cell_node_ids = Table(data,ptr)

  polytopes = fill(QUAD,npanels)
  cell_type = fill(1,npanels)
  reffes = LagrangianRefFE(Float64,QUAD,1)
  cell_reffes=[reffes]

  topo = UnstructuredGridTopology(nodes_3d,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())

  d_to_num_dfaces = [ num_faces(topo,d) for d in 0:num_dims(topo)]
  labels = FaceLabeling(d_to_num_dfaces)

  get_face_entity(labels,0) .= get_isboundary_face(topo,0) #.+ 1
  get_face_entity(labels,1) .= fill(true,length(get_isboundary_face(topo,1))) .+ 1
  get_face_entity(labels,2) .= get_isboundary_face(topo,2) #.+ 1

  add_tag!(labels,"interior",[1,])
  add_tag!(labels,"boundary",[2,])
  add_tag_from_tags!(labels,"all",["interior","boundary"])

  panel_grid = Gridap.Geometry.UnstructuredGrid(nodes_3d,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

  panel_ids = collect(2)
  return panel_grid,topo,labels,panel_ids
end


RADIUS = 1.0*sqrt(3)
a = π/4
cube_grid_3D,topo_3D,face_labels_3D,panel_ids = coarse_panel2_3D(a)
cube_model_3D = UnstructuredDiscreteModel(cube_grid_3D,topo_3D,face_labels_3D)
A_bump, B_bump, b_bump = bump_matrics(a)

writevtk(cube_model_3D,dir*"/panel",append=false)




### apply refinement
cube_model_3D = Gridap.Adaptivity.refine(cube_model_3D)
panel_ids = get_panel_ids(cube_model_3D)
panel_ids = fill(2,4)
writevtk(cube_model_3D,dir*"/panel",append=false)

#### make 2D surface
cube_grid_3D = get_grid(cube_model_3D)

cmaps_3D = get_cell_map(cube_grid_3D)

k = lazy_map(p->  BumpField(A_bump,B_bump,b_bump) ∘ PanelRotationField(rp1_3D[p]), panel_ids)

parametric_cell_map = lazy_map(∘,k,cmaps_3D)

cell_coords_3D = get_cell_coordinates(cube_grid_3D)
cell_coords_3D_p1 = lazy_map(PanelRotationField(rp1_3D[panel_ids[1]]),cell_coords_3D)
cell_coords_2D = lazy_map(BumpField(A_bump,B_bump,b_bump),cell_coords_3D_p1)

nodes_2D = get_nodes_from_coords(cube_grid_3D,cell_coords_2D)

topo_2D = UnstructuredGridTopology(nodes_2D,
    get_cell_node_ids(cube_grid_3D),get_cell_type(cube_grid_3D),[QUAD],Gridap.Geometry.NonOriented())
face_labels_2D = get_face_labeling(cube_model_3D)

cube_grid_2D = Gridap.Geometry.UnstructuredGrid(nodes_2D,
      get_cell_node_ids(cube_grid_3D),get_reffes(cube_grid_3D),get_cell_type(cube_grid_3D),Gridap.Geometry.NonOriented(),
      nothing,parametric_cell_map)

cube_model_2D = UnstructuredDiscreteModel(cube_grid_2D,topo_2D,face_labels_2D)

Ω = Triangulation(cube_model_2D)
Γ = BoundaryTriangulation(cube_model_2D;tags="boundary")
nΓ = get_normal_vector(Γ)
writevtk(Γ,dir*"/panel2",cellfields=["n"=>nΓ],append=false)

writevtk(Ω,dir*"/panel",append=false)


#################################################################################
function coarse_panel12_3D(a::Real)
  npanels = 2

  nodes_3d = a.* [
    Point(1.0, -1.0, -1.0)  # node 1
    Point(1.0, 1.0, -1.0)   # node 2
    Point(1.0, -1.0, 1.0)   # node 3
    Point(1.0, 1.0, 1.0)    # node 4
    Point(-1.0, -1.0, 1.0)  # node 5
    Point(-1.0, 1.0, 1.0)   # node 6
  ]

  ## CCAM panel ordering
  data = [ 1,2,3,4, 3,4,5,6  ] # reorient + rotated

  ptr = generate_ptr(npanels)
  cell_node_ids = Table(data,ptr)

  polytopes = fill(QUAD,npanels)
  cell_type = fill(1,npanels)
  reffes = LagrangianRefFE(Float64,QUAD,1)
  cell_reffes=[reffes]

  topo = UnstructuredGridTopology(nodes_3d,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())

  d_to_num_dfaces = [ num_faces(topo,d) for d in 0:num_dims(topo)]
  labels = FaceLabeling(d_to_num_dfaces)

  get_face_entity(labels,0) .= get_isboundary_face(topo,0) #.+ 1
  get_face_entity(labels,1) .= get_isboundary_face(topo,1) .+ 1
  get_face_entity(labels,2) .= get_isboundary_face(topo,2) #.+ 1

  add_tag!(labels,"interior",[1,])
  add_tag!(labels,"boundary",[2,])
  add_tag_from_tags!(labels,"all",["interior","boundary"])

  panel_grid = Gridap.Geometry.UnstructuredGrid(nodes_3d,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

  panel_ids = collect(1:2)
  return panel_grid,topo,labels,panel_ids
end

RADIUS = 1.0*sqrt(3)
a = π/4
cube_grid_3D,topo_3D,face_labels_3D,panel_ids = coarse_panel12_3D(a)
cube_model_3D = UnstructuredDiscreteModel(cube_grid_3D,topo_3D,face_labels_3D)
A_bump, B_bump, b_bump = bump_matrics(a)

writevtk(cube_model_3D,dir*"/panel",append=false)

cube_model_3D = Gridap.Adaptivity.refine(cube_model_3D)
panel_ids = get_panel_ids(cube_model_3D)
writevtk(cube_model_3D,dir*"/panel",append=false)


cube_grid_3D = get_grid(cube_model_3D)

cmaps_3D = get_cell_map(cube_grid_3D)

k = lazy_map(p-> ShiftField(p)∘  BumpField(A_bump,B_bump,b_bump) ∘ PanelRotationField(rp1_3D[p]), panel_ids)

parametric_cell_map = lazy_map(∘,k,cmaps_3D)

evaluate(parametric_cell_map[2],Point(0,0))

cell_coords_3D = get_cell_coordinates(cube_grid_3D)
# cell_coords_3D_p1 = lazy_map(PanelRotationMap(rp1_3D),cell_coords_3D,panel_ids)
# cell_coords_2D = lazy_map(BumpMap(),cell_coords_3D_p1)

cell_coords_2D = lazy_map(evaluate,k,cell_coords_3D)


nodes_2D = get_nodes_from_coords(cube_grid_3D,cell_coords_2D)

topo_2D = UnstructuredGridTopology(nodes_2D,
    get_cell_node_ids(cube_grid_3D),get_cell_type(cube_grid_3D),[QUAD],Gridap.Geometry.NonOriented())
face_labels_2D = get_face_labeling(cube_model_3D)

cube_grid_2D = Gridap.Geometry.UnstructuredGrid(nodes_2D,
      get_cell_node_ids(cube_grid_3D),get_reffes(cube_grid_3D),get_cell_type(cube_grid_3D),Gridap.Geometry.NonOriented(),
      nothing,parametric_cell_map)

cube_model_2D = UnstructuredDiscreteModel(cube_grid_2D,topo_2D,face_labels_2D)


Ω = Triangulation(cube_model_2D)
Γ = BoundaryTriangulation(cube_model_2D;tags="boundary")
nΓ = get_normal_vector(Γ)
writevtk(Γ,dir*"/panel12",cellfields=["n"=>nΓ],append=false)

writevtk(cube_model_2D,dir*"/panel12",append=false)





##### FE problem
Ω = Triangulation(cube_model_2D)
Γ = BoundaryTriangulation(cube_model_2D;tags="boundary")

_m = map(p->GenericField(met(p)),panel_ids)
M = CellData.GenericCellField(_m,Ω,PhysicalDomain())

_minv = map(p->GenericField(invmet(p)),panel_ids)
invM = CellData.GenericCellField(_minv,Ω,PhysicalDomain())

_msqrt = map(p->GenericField(sqrtmet(p)),panel_ids)
sqrtM = CellData.GenericCellField(_msqrt,Ω,PhysicalDomain())

mNew =  Metric(M,sqrtM,invM,met,sqrtmet,invmet)


dΓ = Measure(Γ,10)
n_Γ = get_normal_vector(Γ)

dΩ = Measure(Ω,10)

ucf = CellField(u_scalar,Ω)
lap_ucf = surface_laplacian(ucf,mNew)
rhs = ucf + lap_ucf

writevtk(Ω,dir*"/panel12_u",cellfields=["u"=>ucf],append=false)


V = TestFESpace(cube_model_2D, ReferenceFE(lagrangian,Float64,1); conformity=:H1)
U = TrialFESpace(V)

_Jcf = map(p->GenericField(Jp(p)),panel_ids)
Jcf = CellData.GenericCellField(_Jcf,Ω,PhysicalDomain())
Gcf = (Operation(transpose)(Jcf)) ⋅ Jcf
invG = Operation(inv)(Gcf)
measG = Operation((meas))(Gcf)
sqrtmeasG = Operation(sqrt)(measG)

mass(u,v) = ∫( (u*v)* (sqrtmeasG) )dΩ
stiffnes(u,v) = ∫( (Jcf ⋅(invG ⋅ ∇(v)))⋅ (Jcf ⋅invG ⋅ ∇(u) ) *(sqrtmeasG))dΩ
bound(v) = ∫( (v*(invG⋅ ∇(ucf))⋅n_Γ )*(sqrtmeasG) )dΓ

helmholtz_biform(u,v) = mass(u,v) - stiffnes(u,v)
helmholtz_liform(v) = ∫(  (rhs*v) * (sqrtmeasG) )dΩ - bound(v)

op = AffineFEOperator(helmholtz_biform,helmholtz_liform,U,V)

uh = solve(LUSolver(),op)

e = l2(uh-ucf,dΩ)

writevtk(Ω,dir*"/panel12_u",cellfields=["u"=>ucf,"uh"=>uh,"e"=>ucf-uh],append=false)
