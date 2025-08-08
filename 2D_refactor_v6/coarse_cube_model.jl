using DrWatson
using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity
using LinearAlgebra

include("swapField.jl")
include("overloads.jl")
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

# panels 1, 2, 3, 4, 5, 6
function coarse_cube_surface_3D(a::Real)
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
  data = [ 1,2,3,4, 5,6,3,4, 7,2,6,4, 7,8,6,5, 1,2,8,7,  1,8,3,5 ]
  # data =  [ 1,2,3,4,  5,6,3,4,  7,8,6,5, 1,2,8,7,  ]

  ptr = generate_ptr(npanels)
  cell_node_ids = Table(data,ptr)

  polytopes = fill(QUAD,npanels)
  cell_type = fill(1,npanels)
  reffes = LagrangianRefFE(Float64,QUAD,1)
  cell_reffes=[reffes]

  topo = UnstructuredGridTopology(nodes_3d,cell_node_ids,cell_type,polytopes,Gridap.Geometry.Oriented())
  labels = FaceLabeling(topo)

  cube_grid = Gridap.Geometry.UnstructuredGrid(nodes_3d,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.Oriented())

  panel_ids = collect(1:6)
  return cube_grid,topo,labels,panel_ids
end

cube_grid,topo,labels,panel_ids = coarse_cube_surface_3D(π/4)
cube_model = UnstructuredDiscreteModel(cube_grid,topo,labels)

cube_model = Gridap.Adaptivity.refine(cube_model)
cube_model = Gridap.Adaptivity.refine(cube_model)
# cube_model = Gridap.Adaptivity.refine(cube_model)
# cube_model = Gridap.Adaptivity.refine(cube_model)

# n = Int(num_cells(cube_model)/4)
# panel_ids = vcat(fill(1,n),fill(2,n),fill(4,n),fill(5,n))
panel_ids = get_panel_ids(cube_model)


## make panel grid
cube_grid = get_grid(cube_model)
cube_topo = get_grid_topology(cube_model)

nodes = get_node_coordinates(cube_grid)
cmaps = get_cell_map(cube_grid)

A1 = [0 1 0
      0 0 1]
A2 = [0 1 0
      1 0 0]
A3 = [1 0 0
      0 0 1]
A4 = [0 -1 0
      0 0 1]
A5 = [0 1 0
      -1 0 0]
A6 = [-1 0 0
      0 0 1]
As = [A1,A2,A3,A4,A5,A6]

swap = lazy_map(p->SwapField(TensorValue(As[p])),panel_ids)

new_cmaps = lazy_map(∘,swap,cmaps)


ref_points = get_cell_ref_coordinates(cube_grid)
lazy_map(evaluate,new_cmaps,ref_points)./1

new_nodes = map(x->Point(x[2],x[3]),nodes)
new_grid = Geometry.UnstructuredGrid(new_nodes,get_cell_node_ids(cube_grid),get_reffes(cube_grid),get_cell_type(cube_grid),OrientationStyle(cube_grid),
                    nothing,new_cmaps)

new_topo = UnstructuredGridTopology(new_nodes,get_cell_node_ids(cube_grid),get_cell_type(cube_topo),get_polytopes(cube_topo),OrientationStyle(cube_topo))
new_labels = FaceLabeling(new_topo)
panel_model = UnstructuredDiscreteModel(new_grid,new_topo,new_labels)



################################################################################
##### Single panel: 1
################################################################################
u(X) = X[1]*X[2]#*X[3]

function cartesian_to_latlon(XYZ)
  X,Y,Z = XYZ
  θ = atan(Y,X)
  ϕ = asin(Z/RADIUS)#atan(Z,sqrt(X^2 + Y^2))
  Point(θ,ϕ)
end
uθϕ(θϕ) = sin(θϕ[2])


function uex(p)
  function _uex(αβ)
    XYZ = forward_map(αβ,p)
    θϕ = cartesian_to_latlon(XYZ)
    uθϕ(θϕ)
  end
end

function panel_meas_metric(p)
  function _panel_meas_metric(αβ)
    meas_metric(αβ,p)
  end
end

function panel_metric(p)
  function _panel_metric(αβ)
    metric(αβ,p)
  end
end

function panel_inv_metric(p)
  function _panel_inv_metric(αβ)
    inv_metric(αβ,p)
  end
end

# function CellData.get_cell_points(trian::Triangulation)
#   cell_ref_coords = get_cell_ref_coordinates(trian)
#   cmaps = get_cell_map(trian)
#   cell_phys_coords = lazy_map(evaluate,cmaps,cell_ref_coords)
#   CellPoint(cell_ref_coords,cell_phys_coords,trian,ReferenceDomain())
# end


p_fe = 2
degree = 6*p_fe

Ω = Triangulation(panel_model)
dΩ = Measure(Ω,degree)
pts = get_cell_points(Ω)
quad_pts = get_cell_points(dΩ)

cell_field = map(p->GenericField(uex(p)),panel_ids)
ucf =  CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

cell_field = map(p->GenericField(panel_meas_metric(p)),panel_ids)
meas_metric_cf =  CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

meas_metric_cf(pts)


cell_field = map(p->GenericField(panel_inv_metric(p)),panel_ids)
inv_metric_cf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

cell_field = map(p->GenericField(panel_metric(p)),panel_ids)
metric_cf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

surflap = 1/meas_metric_cf * divergence( meas_metric_cf *( inv_metric_cf ⋅ gradient(ucf) ) )

surflap(pts)
sum(∫( ucf*meas_metric_cf + surflap*meas_metric_cf  )dΩ)

rhs = ucf + surflap


V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
U = TrialFESpace(V)

l2_biform(u,v) = ∫(u*v*meas_metric_cf)dΩ
l2_liform(v) = ∫(  (ucf*v)*meas_metric_cf )dΩ
op = AffineFEOperator(l2_biform,l2_liform,U,V)
uh_l2 = solve(LUSolver(),op)

poisson_biform(u,v) = ∫(u*v*meas_metric_cf)dΩ -  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_metric_cf )dΩ
# poisson_biform(u,v) = ∫(u*v*meas_metric_cf)dΩ -  ∫((  (J_cf⋅gradient(v))⋅ (J_cf⋅gradient(u)) )*meas_metric_cf )dΩ
poisson_liform(v) = ∫(  (rhs*v)*meas_metric_cf )dΩ
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)
uh = solve(LUSolver(),op)

l2(e,dΩ) = sum(∫( e⋅e )dΩ)
e =  l2(uh-ucf,dΩ)

for pid in collect(1:6)
  mask = panel_ids.==pid
  Ωp = Triangulation(panel_model,mask)
  writevtk(Ωp,dir*"/u$pid",cellfields=["u"=>ucf,"uh"=>uh,"e"=>ucf-uh],append=false)
end

writevtk(Triangulation(cube_model),dir*"/cube_model",append=false)
