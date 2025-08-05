using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity
using LinearAlgebra


function generate_ptr(n)
  nvertices = 4
  ptr  = Vector{Int}(undef,n+1)
  ptr[1]=1
  for i=1:n
    ptr[i+1]=ptr[i]+nvertices
  end
  ptr
end

# only panel 1
function coarse_cube_surface_3D(a::Real)
  npanels = 1

  nodes_3d = a.* [
    Point(1.0, -1.0, -1.0)  # node 1
    Point(1.0, 1.0, -1.0)   # node 2
    Point(1.0, -1.0, 1.0)   # node 3
    Point(1.0, 1.0, 1.0)    # node 4
    # Point(-1.0, -1.0, 1.0)  # node 5
    # Point(-1.0, 1.0, 1.0)   # node 6
    # Point(-1.0, 1.0, -1.0)  # node 7
    # Point(-1.0, -1.0, -1.0) # node 8
  ]

  ## CCAM panel ordering
  data = [ 1,2,3,4]

  ptr = generate_ptr(npanels)
  cell_node_ids = Table(data,ptr)

  polytopes = fill(QUAD,npanels)
  cell_type = fill(1,npanels)
  reffes = LagrangianRefFE(Float64,QUAD,1)
  cell_reffes=[reffes]

  topo = UnstructuredGridTopology(nodes_3d,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
  labels = FaceLabeling(topo)

  cube_grid = Gridap.Geometry.UnstructuredGrid(nodes_3d,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())

  panel_ids = collect(1:6)
  return cube_grid,topo,labels,panel_ids
end

cube_grid,topo,labels,panel_ids = coarse_cube_surface_3D(π/4)
cube_model = UnstructuredDiscreteModel(cube_grid,topo,labels)

cube_model = Gridap.Adaptivity.refine(cube_model)
cube_model = Gridap.Adaptivity.refine(cube_model)
cube_model = Gridap.Adaptivity.refine(cube_model)
cube_model = Gridap.Adaptivity.refine(cube_model)

## make panel grid
cube_grid = get_grid(cube_model)
cube_topo = get_grid_topology(cube_model)

nodes = get_node_coordinates(cube_grid)
cmaps = get_cell_map(cube_grid)

A1 = [0 1 0
    0 0 1]

swap = fill(SwapField(TensorValue(A1)),length(cmaps))
new_cmaps = lazy_map(∘,swap,cmaps)

evaluate(new_cmaps[1],Point(1,1))
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
u(X) = X[1]*X[2]*X[3]
## force the panel number below
function uex(αβ)
  XYZ = forward_map(αβ,1)
  u(XYZ)
end

function panel1_meas_metric(αβ)
  meas_metric(αβ,1)
end

function panel1_inv_metric(αβ)
  inv_metric(αβ,1)
end

p = 2
degree = 4*p

Ω = Triangulation(panel_model)
dΩ = Measure(Ω,degree)
pts = get_cell_points(Ω)

ucf = CellField(uex,Ω)
meas_metric_cf = CellField(panel1_meas_metric,Ω)
inv_metric_cf = CellField(panel1_inv_metric,Ω)

gradient(ucf)(pts)

surflap = 1/meas_metric_cf * divergence( meas_metric_cf *( inv_metric_cf ⋅ gradient(ucf) ) )

sum(∫( ucf*meas_metric_cf + surflap*meas_metric_cf  )dΩ)

rhs = ucf + surflap


V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p); conformity=:H1)
U = TrialFESpace(V)

poisson_biform(u,v) = ∫(u*v*meas_metric_cf)dΩ -  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_metric_cf )dΩ
poisson_liform(v) = ∫(  (rhs*v)*meas_metric_cf )dΩ
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

uh = solve(LUSolver(),op)

l2(e,dΩ) = sum(∫( e⋅e )dΩ)
e =  l2(uh-ucf,dΩ)
