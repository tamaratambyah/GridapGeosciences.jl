using DrWatson
using Gridap
using Gridap.Geometry, Gridap.Fields, Gridap.Arrays, Gridap.CellData, Gridap.ReferenceFEs
using Gridap.Adaptivity
using LinearAlgebra

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

# panels 1, 3
function coarse_cube_surface_3D(a::Real)
  npanels = 2

  nodes_3d = a.* [
    Point(1.0, -1.0, -1.0)  # node 1
    Point(1.0, 1.0, -1.0)   # node 2
    Point(1.0, -1.0, 1.0)   # node 3
    Point(1.0, 1.0, 1.0)    # node 4
    Point(-1.0, -1.0, 1.0)  # node 5
    # Point(-1.0, 1.0, 1.0)   # node 6
    # Point(-1.0, 1.0, -1.0)  # node 7
    Point(-1.0, -1.0, -1.0) # node 8
  ]

  ## CCAM panel ordering
  data = [ 1,2,3,4, 1,3,6,5 ]

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
# cube_model = Gridap.Adaptivity.refine(cube_model)
# cube_model = Gridap.Adaptivity.refine(cube_model)

## make panel grid
panel_ids = get_panel_ids(cube_model)
cube_grid = get_grid(cube_model)
cube_topo = get_grid_topology(cube_model)

nodes = get_node_coordinates(cube_grid)
cmaps = get_cell_map(cube_grid)

A1 = [0 1 0
      0 0 1]
A6 = [0 0 1
      -1 0 0]
As = [A1,A6]

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

n = Int(num_cells(cube_model)/2)
panel_ids = vcat(fill(1,n),fill(6,n))



################################################################################
##### Panels 1, 6
################################################################################

function cartesian_to_latlon(XYZ)
  X,Y,Z = XYZ
  θ = rem2pi(atan(Y,X), RoundNearest)
  ϕ = atan(Z,sqrt(X^2 + Y^2))
  Point(θ,ϕ)
end
uθϕ(θϕ) = sin(θϕ[2])

## force the panel number below
function uex(p)
  function _uex(αβ)
    # XYZ = forward_map(αβ,p)
    # θϕ = cartesian_to_latlon(XYZ)
    # uθϕ(θϕ)
    αβ[1]*αβ[2]
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

function CellData.get_cell_points(trian::Triangulation)
  cell_ref_coords = get_cell_ref_coordinates(trian)
  cmaps = get_cell_map(trian)
  cell_phys_coords = lazy_map(evaluate,cmaps,cell_ref_coords)
  CellPoint(cell_ref_coords,cell_phys_coords,trian,ReferenceDomain())
end


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


cell_field = map(p->GenericField(panel_inv_metric(p)),panel_ids)
inv_metric_cf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

cell_field = map(p->GenericField(panel_metric(p)),panel_ids)
metric_cf = CellData.GenericCellField(cell_field,Ω,PhysicalDomain())

surflap = 1/meas_metric_cf * divergence( meas_metric_cf *( inv_metric_cf ⋅ gradient(ucf) ) )

sum(∫( ucf*meas_metric_cf + surflap*meas_metric_cf  )dΩ)

rhs = ucf + surflap

V = TestFESpace(Ω, ReferenceFE(lagrangian,Float64,p_fe); conformity=:H1)
U = TrialFESpace(V)

poisson_biform(u,v) = ∫(u*v*meas_metric_cf)dΩ -  ∫( ( gradient(v)⋅ (inv_metric_cf⋅ gradient(u) ) )*meas_metric_cf )dΩ
poisson_liform(v) = ∫(  (rhs*v)*meas_metric_cf )dΩ
op = AffineFEOperator(poisson_biform,poisson_liform,U,V)

uh = solve(LUSolver(),op)

l2(e,dΩ) = sum(∫( e⋅e )dΩ)
e =  l2(uh-ucf,dΩ)



inv_mapping =  map(p-> InverseMapField2(p) , panel_ids)

cf_mapped = lazy_map(Broadcasting(∘),get_data(ucf),inv_mapping)
cf_ambient = CellData.GenericCellField(cf_mapped,Triangulation(cube_model),PhysicalDomain() )

cf_ambient(get_cell_points(Triangulation(cube_model)))

writevtk(Triangulation(cube_model),dir*"/cube_model",cellfields=["u"=>cf_ambient],append=false)

function inverse_map(p)
  function _inverse_map(XYZ)
    X,Y,Z = XYZ
    α,β = -10.0,-10.0

    if p == 1
      α = atan(Y,X)
      β = atan(Z,X)
    elseif p == 2
      α = atan(Y,Z)
      β = atan(X,Z)
    elseif p == 3
      α = atan(X,Y)
      β = atan(Z,Y)
    elseif p == 4
      α = atan(-Z,X)
      β = atan(Y,X)
    elseif p == 5
      α = atan(X,-Z)
      β = atan(Y,Z)
    elseif p == 6
      α = atan(-Z,Y)
      β = atan(X,Y)
    end

    αβ = Point(α,β)
    return αβ
  end

end

for pid in collect(1:6)
  mask = panel_ids.==pid
  Ωp = Triangulation(panel_model,mask)
  writevtk(Ωp,dir*"/u$pid",cellfields=["u"=>ucf,"uh"=>uh,"e"=>ucf-uh],append=false)
end

writevtk(Triangulation(cube_model),dir*"/cube_model",append=false)
