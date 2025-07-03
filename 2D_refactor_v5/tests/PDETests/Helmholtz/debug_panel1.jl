"""
build panel 1 by hand
"""

using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../../../src/initialise.jl")

function coarse_panel_3D(a::Real)
  npanels = 1

  nodes_3d = a.* [
    Point(1.0, -1.0, -1.0)  # node 1
    Point(1.0, 1.0, -1.0)   # node 2
    Point(1.0, -1.0, 1.0)   # node 3
    Point(1.0, 1.0, 1.0)    # node 4
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
  get_face_entity(labels,1) .= get_isboundary_face(topo,1) .+ 1
  get_face_entity(labels,2) .= get_isboundary_face(topo,2) #.+ 1

  add_tag!(labels,"interior",[1,])
  add_tag!(labels,"boundary",[2,])
  add_tag_from_tags!(labels,"all",["interior","boundary"])

  panel_grid = Gridap.Geometry.UnstructuredGrid(nodes_3d,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

  panel_ids = collect(1)
  return panel_grid,topo,labels,panel_ids
end


RADIUS = 1.0*sqrt(3)
a = π/4
cube_grid_3D,topo_3D,face_labels_3D,panel_ids = coarse_panel_3D(a)
cube_model_3D = UnstructuredDiscreteModel(cube_grid_3D,topo_3D,face_labels_3D)
A_bump, B_bump, b_bump = bump_matrics(a)


### apply refinement
cube_model_3D = Gridap.Adaptivity.refine(cube_model_3D)
panel_ids = get_panel_ids(cube_model_3D)

#### make 2D surface
cube_grid_3D = get_grid(cube_model_3D)
cmaps_3D = get_cell_map(cube_grid_3D)

k = lazy_map(p->  BumpField(A_bump,B_bump,b_bump) ∘ PanelRotationField(rp1_3D[p]), panel_ids)

parametric_cell_map = lazy_map(∘,k,cmaps_3D)

cell_coords_3D = get_cell_coordinates(cube_grid_3D)
cell_coords_2D = lazy_map(BumpField(A_bump,B_bump,b_bump),cell_coords_3D)

nodes_2D = get_nodes_from_coords(cube_grid_3D,cell_coords_2D)

topo_2D = UnstructuredGridTopology(nodes_2D,
    get_cell_node_ids(cube_grid_3D),get_cell_type(cube_grid_3D),[QUAD],Gridap.Geometry.NonOriented())
face_labels_2D = get_face_labeling(cube_model_3D)

cube_grid_2D = Gridap.Geometry.UnstructuredGrid(nodes_2D,
      get_cell_node_ids(cube_grid_3D),get_reffes(cube_grid_3D),get_cell_type(cube_grid_3D),Gridap.Geometry.NonOriented(),
      nothing,parametric_cell_map)

cube_model_2D = UnstructuredDiscreteModel(cube_grid_2D,topo_2D,face_labels_2D)
writevtk(cube_model_2D,dir*"/panel",append=false)


function uθϕ(αβ)
    θϕ = GnomonicField()(αβ)
    cos(θϕ[1])*sin(θϕ[2])
end


##### FE problem
Ω = Triangulation(cube_model_2D)
Γ = BoundaryTriangulation(cube_model_2D;tags="boundary")

m = Metric(cubedsphere,Ω)
mΓ = Metric(cubedsphere,Γ)

dΓ = Measure(Γ,10)
nΓ = get_normal_vector(Γ)

dΩ = Measure(Ω,10)

ucf = CellField(uθϕ,Ω)
lap_ucf = surface_laplacian(ucf,m)
rhs = ucf + lap_ucf
sum(∫(rhs)dΩ)

V = TestFESpace(cube_model_2D, ReferenceFE(lagrangian,Float64,2); conformity=:H1)
U = TrialFESpace(V)

M = CellField(metric1,Ω)
invM = CellField(invmetric1,Ω)
sqrtM = CellField(sqrtmet1,Ω)

mass(u,v) = ∫( (u*v)*sqrtM )dΩ
stiffnes(u,v) = ∫( (∇(v)⋅ (invM⋅ ∇(u) )) *sqrtM)dΩ
bound(v) = ∫( (((invM ⋅∇(ucf))⋅(nΓ)*v)*sqrtM ) )dΓ
force(v) = ∫(  (rhs*v)*sqrtM )dΩ

bb = assemble_vector(bound,V)

helmholtz_biform(u,v) = mass(u,v) - stiffnes(u,v)
helmholtz_liform(v) = force(v) #- bound(v)

op = AffineFEOperator(helmholtz_biform,helmholtz_liform,U,V)

uh = solve(LUSolver(),op)

e = l2(uh-ucf,dΩ)

writevtk(Ω,dir*"/poisson_panel",cellfields=["u"=>ucf,"uh"=>uh,"e"=>ucf-uh],append=false)
