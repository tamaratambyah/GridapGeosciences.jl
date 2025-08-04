using DrWatson
using Gridap
using Gridap.Arrays, Gridap.ReferenceFEs, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.Adaptivity, Gridap.Helpers, Gridap.TensorValues, Gridap.Fields
using Test, BenchmarkTools
using LinearAlgebra
using FillArrays

include("../../../src/initialise.jl")

pt1 = Point(π/8,π/4)
Jpanel1(pt1)

pt2 = Point(π/8,-π/4)
Jpanel2(pt2)

RADIUS = sqrt(3)


P1 = one(TensorValue{2,2,Float64})
P2 = TensorValue{2,2,Float64}(1,0.0,0.0,-1)
P3 = TensorValue{2,2,Float64}(-1,0,0,1)
P4 = TensorValue{2,2,Float64}(0,1,-1,0)
P5 = TensorValue{2,2,Float64}(0,1,1,0)
P6 = TensorValue{2,2,Float64}(0,1,1,0)

Rs_p1 = [rp1_3D[1];rp1_3D[2];rp1_3D[3];rp1_3D[4];rp1_3D[5];rp1_3D[6]]
Rs_1p = map(x->inv(x),Rs_p1) #[r1p_3D[1];r1p_3D[2];r1p_3D[3];r1p_3D[4];r1p_3D[5];r1p_3D[6]]
Ps = [P1;P2;P3;P4;P5;P6]
Psinv = map(x->inv(x),Ps) #[P1;inv(P2);inv(P3);inv(P5);inv(P6)]



function coarse_panels_13_3D(a::Real)
  npanels = 6

  topo_nodes = a.* [
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
  vertex_data = [ 1,2,3,4, 3,4,5,6, 2,7,4,6, 8,5,7,6, 1,8,2,7, 1,3,8,5  ]
  # vertex_data = [ 1,2,3,4, 3,4,5,6, 2,7,4,6,  1,8,2,7, 1,3,8,5  ]

  ptr = generate_ptr(npanels)
  cell_vertices = Table(vertex_data,ptr)

  polytopes = fill(QUAD,npanels)
  cell_type = fill(1,npanels)
  reffes = LagrangianRefFE(Float64,QUAD,1)
  cell_reffes=[reffes]

  topo = UnstructuredGridTopology(topo_nodes,cell_vertices,cell_type,polytopes,Gridap.Geometry.NonOriented())

  d_to_num_dfaces = [ num_faces(topo,d) for d in 0:num_dims(topo)]
  labels = FaceLabeling(d_to_num_dfaces)

  get_face_entity(labels,0) .= get_isboundary_face(topo,0) #.+ 1
  get_face_entity(labels,1) .= get_isboundary_face(topo,1) .+ 1
  get_face_entity(labels,2) .= get_isboundary_face(topo,2) #.+ 1

  add_tag!(labels,"interior",[1,])
  add_tag!(labels,"boundary",[2,])
  add_tag_from_tags!(labels,"all",["interior","boundary"])


  grid_nodes = a.* [
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
  grid_data = [ 1,2,3,4, 3,4,5,6, 2,7,4,6, 8,5,7,6, 1,8,2,7, 1,3,8,5  ]
  # grid_data = [ 1,2,3,4, 3,4,5,6, 2,7,4,6, 1,8,2,7, 1,3,8,5  ]

  ptr = generate_ptr(npanels)
  cell_node_ids = Table(grid_data,ptr)

  panel_grid = Gridap.Geometry.UnstructuredGrid(grid_nodes,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

  panel_ids = collect(1:3)
  return panel_grid,topo,labels,panel_ids
end

A_bump, B_bump, b_bump =  bump_matrics(π/4)

panel_grid,topo,labels,panel_ids = coarse_panels_13_3D(π/4)

_model = UnstructuredDiscreteModel(panel_grid,topo,labels)

_model = Gridap.Adaptivity.refine(_model)
writevtk(_model,dir*"/panel13",append=false)

panel_ids = get_panel_ids(_model)

### make parametric grid
grid = get_grid(_model)
topo = get_grid_topology(_model)
cmaps = get_cell_map(grid)


k = lazy_map(p->  InversionField(Ps[p]) ∘ BumpField(A_bump,B_bump,b_bump) ∘ PanelRotationField(Rs_p1[p]), panel_ids)
_cmaps = lazy_map(∘,k,cmaps)

cell_coords = lazy_map(evaluate,_cmaps,get_cell_ref_coordinates(grid))
nodes = get_panel_1_nodes_from_coords(grid,cell_coords,panel_ids)

_grid = UnstructuredGrid(nodes,get_cell_node_ids(grid),get_reffes(grid),get_cell_type(grid),
  OrientationStyle(grid),nothing,_cmaps)

nodes_topo = get_panel_1_nodes_from_coords(topo,cell_coords,panel_ids)
_topo = UnstructuredGridTopology(nodes,get_faces(topo,2,0),get_cell_type(topo),get_polytopes(topo),Gridap.Geometry.NonOriented())

model = UnstructuredDiscreteModel(_grid,_topo,FaceLabeling(_topo))

writevtk(model,dir*"/panel13",append=false)

######### FE
Ω = Triangulation(model)
writevtk(Ω,dir*"/test",append=false)

dΩ = Measure(Ω,10)


V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,3); conformity=:H1)
U = TrialFESpace(V)

uX_scalar(x) = x[1]*x[2]*x[3]
uθϕ_scalar(x) =  cos(x[1])*sin(x[2])
function u_scalar_ambient2parametric2(p::Int,uX::Function)
  function _u(αβ)
    _αβ =  InversionField(Psinv[p])(αβ)
    θϕ = GnomonicField()(_αβ)
    _XYZ = SigmaField(RADIUS)(θϕ)
    XYZ = PanelRotationField(Rs_1p[p])(_XYZ)
    uX(XYZ)
    # θϕp = InvSigmaField(RADIUS)(XYZ)
    # uX(θϕp)
  end
end

cell_field = map(p->GenericField(u_scalar_ambient2parametric2(p,uX_scalar)),panel_ids)
# cell_field = map(p->GenericField(u_scalar_ambient2parametric2(p,uθϕ_scalar)),panel_ids[panel_ids.==1])
ucf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())
writevtk(Ω,dir*"/test_u",cellfields=["u"=>ucf],append=false)


M = CellField(metric1,Ω)
invM = CellField(invmetric1,Ω)
sqrtM = CellField(sqrtmet1,Ω)

mNew =  Metric(M,sqrtM,invM,metric1,sqrtmet1,invmetric1)

lapucf = surface_laplacian(ucf,mNew) #1/sqrtM * divergence( sqrtM*( invM ⋅ gradient(ucf) ) )

rhs = ucf + lapucf
sum(∫(rhs)dΩ)

mass(u,v) = ∫( (u*v)*sqrtM )dΩ
stiffnes(u,v) = ∫( (∇(v)⋅ (invM⋅ ∇(u) )) *sqrtM)dΩ
# bound(v) = ∫( (((invM ⋅∇(ucf))⋅(nΓ)*v)*sqrtM ) )dΓ
force(v) = ∫(  (rhs*v)*sqrtM )dΩ

helmholtz_biform(u,v) = mass(u,v) - stiffnes(u,v)
helmholtz_liform(v) = force(v) #- bound(v)

op = AffineFEOperator(helmholtz_biform,helmholtz_liform,U,V)

A = assemble_matrix(mass,U,V)
println(diag(Array(A)))
b = assemble_vector(force,V)

uh = solve(LUSolver(),op)

e = l2(uh-ucf,dΩ)

# Ωp = Triangulation(model,panel_ids.==4)
writevtk(Ω,dir*"/test_u",cellfields=["u"=>ucf,"uh"=>uh,"e"=>ucf-uh,],append=false)

######## map to proper parametric space

inv_mapping = lazy_map(p->  InversionField(Ps[p]) ∘ BumpField(A_bump,B_bump,b_bump) ∘ PanelRotationField(Rs_p1[p]), panel_ids)
_Ω = Triangulation(_model)
_pts = get_cell_points(_Ω)

_cf  = change_domain(uh.cell_field,ReferenceDomain(),PhysicalDomain())
cf_mapped = lazy_map(Broadcasting(∘),get_data(_cf),inv_mapping)
uh_mapped = CellData.GenericCellField(cf_mapped,_Ω,PhysicalDomain() )
# uh_mapped(_pts)

_cf_mapped = lazy_map(Broadcasting(∘),get_data(ucf),inv_mapping)
ucf_mapped = CellData.GenericCellField(_cf_mapped,_Ω,PhysicalDomain() )

# ucf_mapped(_pts)

_Ωp = Triangulation(_model,panel_ids.!=4)
writevtk(_Ωp,dir*"/test_u_mapped",
cellfields=["u"=>ucf_mapped,"uh"=>uh_mapped,"e"=>ucf_mapped-uh_mapped],append=false)

# writevtk(_Ω,dir*"/test_u_mapped",
# cellfields=["u"=>ucf_mapped],append=false)
