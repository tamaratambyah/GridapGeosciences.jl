using Gridap
using GridapGeosciences
using Gridap.Geometry, Gridap.Arrays, Gridap.ReferenceFEs
using Gridap.Adaptivity, Gridap.Fields
import GridapGeosciences.Geometry: _CCAM_cube_nodes_3d
import GridapGeosciences.Helpers: normal_vec

using DrWatson
dir = datadir("Extruded_model")
!isdir(dir) && mkdir(dir)

extrusion(x) = x +  sqrt(3)*normal_vec(x)
_intrusion(x) = x - sqrt(3)*normal_vec(x)

a = π/4
x = a*Point(1.0,-1.0,-1.0)
ex = extrusion(x)
_intrusion(ex)
axis = Point(0.0,-1.0,0.0)
proj = (ex⋅axis)
γ = proj - a



function generate_ptr(n)
  nvertices = 8
  ptr  = Vector{Int}(undef,n+1)
  ptr[1]=1
  for i=1:n
    ptr[i+1]=ptr[i]+nvertices
  end
  ptr
end

npanels = 6

data = [ 1,2,3,4,  9,10,11,12,
         3,4,5,6,  11,12,13,14,
         2,7,4,6,  10,15,12,14,
         8,5,7,6,  16,13,15,14,
         1,8,2,7,  9,16,10,15,
         1,3,8,5,  9,11,16,13  ]
ptr = generate_ptr(npanels)
cell_node_ids = Table(data,ptr)

polytopes = fill(HEX,1)
cell_type = fill(1,npanels)
reffes = LagrangianRefFE(Float64,HEX,1)
cell_reffes=[reffes]

a = π/4 ### for testing, use a = 1.0
cube_surface_nodes = _CCAM_cube_nodes_3d(a)

extrusion(γ::Float64,x) = x + γ*sqrt(3)*normal_vec(x)
_intrusion(γ::Float64,x) = x - γ*sqrt(3)*normal_vec(x)

extrusion_nodes_3D =  [extrusion.(0.0,cube_surface_nodes)
                      extrusion.(1.0,cube_surface_nodes)]

topo = UnstructuredGridTopology(extrusion_nodes_3D,cell_node_ids,cell_type,polytopes,Gridap.Geometry.NonOriented())
labels = FaceLabeling(topo)
grid = Gridap.Geometry.UnstructuredGrid(extrusion_nodes_3D,cell_node_ids,cell_reffes,cell_type,Gridap.Geometry.NonOriented())
extruded_cube_model = UnstructuredDiscreteModel(grid,topo,labels)
extruded_cube_model_ref = Adaptivity.refine(extruded_cube_model)

writevtk(Triangulation(extruded_cube_model),dir*"/extruded_cube",append=false)
writevtk(Triangulation(extruded_cube_model_ref),dir*"/extruded_cube_ref",append=false)

panel_normal = [
  Point(1.0, 0.0, 0.0) #+X
  Point(0.0, 0.0, 1.0) #+Z
  Point(0.0, 1.0, 0.0) #+Y
  Point(-1.0, 0.0, 0.0) #-X
  Point(0.0, 0.0, -1.0) #-Z
  Point(0.0, -1.0, 0.0) #-Y
]


function extrusion_variable(p::Int,x::VectorValue{3};a=π/4)
  n = panel_normal[p]
  proj = (x⋅n)
  γ = proj - a
  return γ
end

x = extrusion_nodes_3D[1:4]
γ =  extrusion_variable.(1,x)
_intrusion.(γ,x)

function map_points_to_surface(p::Int,x::VectorValue{3})
  γ = extrusion_variable(p,x)
  _intrusion(γ,x)
end


As = map(p->extruded_cube_to_αβγ(p),collect(1:6))
b = VectorValue(0.0,0.0,1.0)
p = 1
f = Map2ExtrudedPanel(p,As[p],b)

x = extrusion_nodes_3D[9]
f(x)


struct Map2ExtrudedPanel{A,B,C}  <: Field
  p::A
  Apanel::B
  b::C # z vector (0,0,1)
end


function Gridap.Arrays.return_cache(f::Map2ExtrudedPanel,cellx::AbstractArray{<:VectorValue{3}})
  x = first(cellx)
  T = typeof(x)
  y = similar(cellx,T)
  γ = similar(cellx,Float64)
  return y,γ
end


function Gridap.Arrays.evaluate!(cache,f::Map2ExtrudedPanel,cellx::AbstractArray{<:VectorValue{3}} )
  p = f.p
  y,γ = cache
  map!(x -> map_points_to_surface(p,x), y, cellx)
  return y
end


function Gridap.Arrays.return_cache(f::Map2ExtrudedPanel,x::VectorValue{3})
  T = typeof(x)
  y = zero(T)
  γ = 0.0
  return y,γ
end

function Gridap.Arrays.evaluate!(cache,f::Map2ExtrudedPanel,x::VectorValue{3})
  p = f.p
  Apanel = f.Apanel
  b = f.b
  y,γ = cache
  γ = extrusion_variable(p,x)
  y = _intrusion(γ,x) # surface points
  return Apanel ⋅ y + γ*b
end


function extruded_cube_to_αβγ(p::Int)
  if p == 1
    A = [0 1 0
        0 0 1
        0 0 0] # γ is X component
  elseif p == 2
    A = [0 1 0
        -1 0 0
        0 0 0] # γ is Z component
  elseif p == 3
    A = [-1 0 0
          0 0 1
          0 0 0]  # γ is Y component
  elseif p == 4
    A = [0 0 1
         0 1 0
         0 0 0] # γ is X component
  elseif p == 5
    A = [-1 0 0
          0 1 0
         0 0 0]  # γ is Z component
  elseif p == 6
    A = [0 0 1
        -1 0 0
         0 0 0]  # γ is Y component
  end

  return TensorValue(A)

end


model = extruded_cube_model

panel_ids = get_panel_ids(model)
extruded_cube_cmaps = get_cell_map(get_grid(model))
A_extruedcube_2_extrudedpanel = map(p->extruded_cube_to_αβγ(p),collect(1:6))


## map all 3D points-> surface
##    * in this step, store the extrusion variable
## map surface points -> alpha,beta
## add back the extrusion variable
extrudedcube2panel_map = lazy_map(p->MatMultField( A_extruedcube_2_extrudedpanel[p] ) ∘ Map2Surface(p), panel_ids)

extruded_panel_cmaps = lazy_map(∘,extrudedcube2panel_map,extruded_cube_cmaps)
ref_pts = get_cell_ref_coordinates(get_grid(model))
extruded_panel_coords = lazy_map(evaluate,extruded_panel_cmaps,ref_pts)


####### extruded panel model
extruded_cube_grid = get_grid(model)
extruded_topo = get_grid_topology(model)

extruded_panel_nodes = get_node_coordinates(extruded_cube_grid) # these are just junk nodes, never used
extruded_panel_grid = Geometry.UnstructuredGrid(extruded_panel_nodes,
    get_cell_node_ids(extruded_cube_grid),get_reffes(extruded_cube_grid),
    get_cell_type(extruded_cube_grid),OrientationStyle(extruded_cube_grid),
    nothing,extruded_panel_cmaps)
extruded_panel_topo = UnstructuredGridTopology(extruded_panel_nodes,get_cell_node_ids(extruded_cube_grid),
      get_cell_type(extruded_topo),get_polytopes(extruded_topo),OrientationStyle(extruded_topo))
extruded_panel_labels = FaceLabeling(extruded_panel_topo)

extruded_panel_model = ParametricDiscreteModel(extruded_panel_grid,extruded_panel_topo,extruded_panel_labels,panel_ids)

include("forward_map_3D.jl")
cell_geo_map = geo_map_func_3D(get_panel_ids(extruded_panel_model))
writevtk(Triangulation(extruded_panel_model),dir*"/extruded_panel_model",append=false,geo_map=cell_geo_map)

extruded_model = lazy_map(∘,cell_geo_map,extruded_panel_cmaps)
lazy_map(evaluate,extruded_model,ref_pts)

p = 6
pts = extruded_panel_coords[p][5]
α,β,γ = pts

#### compute XYZ point on surface of inner sphere using 2D forward_map
αβ = Point(α,β)
XYZ_surf = forward_map(p,αβ)
XYZ_surf + normal_vec(XYZ_surf)
